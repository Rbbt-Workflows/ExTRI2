'''Run from postprocessing_checkings.ipynb'''

import pandas as pd
import numpy as np
import os
import json
from Bio import Entrez
import time
import difflib
import re
import requests
import time

from scripts.postprocessing.renormalisations import renormalize, discard, fix_NFKB_AP1

# *Always* tell NCBI who you are
Entrez.email = "example24@gmail.com"

TaxID = {
    'human': '9606',
    'mouse': '10090',
    'rat':   '10116',
}

def load_config() -> dict:
    config = {
        "data_p": "../../data/",
        "temp_p": "../../data/tmp/",
        "results_p": "../../results/",
        "chunk_size": 1_500_000
    }
    RESULTS_P = config['results_p']
    DATA_P = config['data_p']

    # Inputs
    config['raw_results_p'] = RESULTS_P + 'ExTRI2.tsv'
    config['orthologs_p'] = DATA_P + 'external/human_HGNC_orthologs/'
    
    # External
    config['TF_types'] = ['dbTF', 'GTF', 'coTF', 'll_coTF']
    config['TF_ids_p'] = DATA_P + "postprocessing/"
    config['ensembl_folder'] = DATA_P + 'external/ensembl_Release_115_orthologs/'

    # Temporary
    TEMP_P = config['temp_p']
    config['raw_TRI_p']                   = TEMP_P + 'ExTRI2_valid_raw.tsv'
    config['raw_nonTRI_sample_p']         = TEMP_P + 'ExTRI2_nonvalid_subset_raw.tsv'
    config['TRI_pre_renorm_p']            = TEMP_P + 'ExTRI2_valid_prenorm.tsv'               # Used for validation purposes
    config['nonTRI_sample_p']             = TEMP_P + 'ExTRI2_nonvalid_subset.tsv'             # Used for validation purposes
    config['EntrezID_to_Symbol_TRI_p']    = TEMP_P + 'EntrezID_to_Symbol_TaxID.json'
    config['EntrezID_to_Symbol_nonTRI_p'] = TEMP_P + 'EntrezID_to_Symbol_TaxID_nonvalid.json'

    # Tables
    POSTPROCESS_TABLES_P = DATA_P + 'postprocessing/tables/'
    config['postprocess_tables_folder']  = POSTPROCESS_TABLES_P
    config['discarded_sents_p']          = POSTPROCESS_TABLES_P + 'discarded_sents.tsv'
    config['renormalized_sents_p']       = POSTPROCESS_TABLES_P + 'renormalized_sents.tsv'
    config['orthologs_final_p']          = POSTPROCESS_TABLES_P + 'orthologs_final.tsv'   
    
    config['NFKB_AP1_discarded_sents_p'] = POSTPROCESS_TABLES_P + 'NFKB_AP1_discarded_sents.tsv'
    config['AP1_NFKB_breakdown_p']       = POSTPROCESS_TABLES_P + 'AP1_NFKB_breakdown.tsv'
    config['AP1_NFKB_breakdown_cols']    = ['entity', 'dimer', 'symbol', '% unmodified', '% modified', 'regex unmodified', 'regex modified']

    # Outputs
    config['final_ExTRI2_p'] = RESULTS_P + "ExTRI2_final_resource.tsv"

    return config

def load_df(df_path: str, header=0) -> pd.DataFrame:
    df = pd.read_csv(df_path, sep='\t', header=header, keep_default_na=False, dtype='str')
    return df

def retrieve_Entrez_annotations(id_list: list):
    """
    Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    
    Returns the dictionary: {ID: {'Name': "", 'TaxID': ""}}
    """

    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)
    webEnv = result["WebEnv"] # type: ignore
    queryKey = result["QueryKey"] # type: ignore
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
    annotationsSummary = annotations['DocumentSummarySet']['DocumentSummary'] # type: ignore

    assert len(id_list) == len(annotationsSummary), f"id_list and annotationsSummary are of different length: {len(id_list)} != {len(annotationsSummary)}"

    return annotationsSummary

def get_all_ids(df) -> set:
    '''Get list of unique NCBI IDs for TF and TG''' 
    TF_IDs = set(id for ids in df['TF Id'] for id in ids.split(";"))
    TG_IDs = set(id for ids in df['TG Id'] for id in ids.split(";"))
    IDs = (TF_IDs | TG_IDs)
    IDs -= {'', 'nan', 'None'}
    assert all(id.isdigit() for id in IDs), 'Some IDs are not digits'
    print(f"We got {len(TF_IDs)} different TFs and {len(TG_IDs)} different TGs from sentences labelled as TRI")

    return IDs

def save_Symbol_TaxID_dict(df: pd.DataFrame, output_path: str) -> None:
    '''
    Given a dataframe with the columns "TF Id" and "TG Id", retrieve its symbol and TaxID from Entrez and 
    save a json dictionary with the structure:
    { 
    <EntrezID>: { Name:  <Symbol>, TaxID: <TaxID> }
    ...
    }
    '''

    IDs = get_all_ids(df)

    # If dict already exists, eliminate already downloaded IDs (not to download them again)
    if os.path.exists(output_path):

        with open(output_path, 'r') as f:
            past_EntrezID_to_Symbol = json.load(f)
        ids_in_dict = past_EntrezID_to_Symbol.keys()
        IDs -= ids_in_dict
    else:
        past_EntrezID_to_Symbol = {}
    
    # Retrieve from Entrez in blocks of 9000
    print('Retrieving from Entrez...\n')
    annotationsSummary = []
    IDs_list = list(IDs)
    total_IDs = len(IDs_list)

    for i in range(0, total_IDs, 9000):
        IDs_subset = IDs_list[i:i+9000]
        annotationsSummary_subset = retrieve_Entrez_annotations(IDs_subset)
        annotationsSummary.extend(annotationsSummary_subset)
        print(f'Retrieved {len(annotationsSummary)}/{total_IDs} IDs from Entrez', end='\r')
        time.sleep(5)

    # Get relevant info
    EntrezID_to_Symbol = {ID : {'Name': annotation["Name"], 'TaxID': annotation['Organism']['TaxID']} for ID, annotation in zip(IDs, annotationsSummary)}
    EntrezID_to_Symbol[""] = {'Name': '', 'TaxID': ''}
    EntrezID_to_Symbol['nan'] = {'Name': '', 'TaxID': ''}

    # Join w/ previously downloaded
    EntrezID_to_Symbol = {**EntrezID_to_Symbol, **past_EntrezID_to_Symbol}

    # Save
    with open(output_path, 'w') as f:
        json.dump(EntrezID_to_Symbol, f)

    return 

def get_TRI_nonTRI_raw_dbs(ExTRI2_results_path: str, raw_TRI_path: str, raw_nonTRI_sample_path: str, chunk_size = 1_200_000) -> None: 
    '''
    ExTRI2 contains millions of candidate sentences.
    This function creates 2 more manageable datasets: one with all sentences labelled as TRI, one with a sample of non-TRI sentences
    '''
    TRI_df_list = []
    nonTRI_df_sample_list = []

    # Load the sentences in chunks of chunk_size
    candidate_sents = pd.read_csv(ExTRI2_results_path, sep='\t', header=1, chunksize=chunk_size, keep_default_na=False)

    for i, chunk in enumerate(candidate_sents):

        # "TRI" is labelled as "Valid" in the raw ExTRI2 file
        TRI_df = chunk[chunk['Valid'] == 'Valid']
        nonTRI_df_sample = chunk[chunk['Valid'] == 'Non valid'].sample(n=5000)

        TRI_df = TRI_df.rename(columns={'Text': 'Sentence', 'Gene': 'TG', 'IG Id': 'TG Id'})
        nonTRI_df_sample = nonTRI_df_sample.rename(columns={'Text': 'Sentence', 'Gene': 'TG', 'IG Id': 'TG Id'})

        TRI_df_list.append(TRI_df)
        nonTRI_df_sample_list.append(nonTRI_df_sample)

        print(f'Processed {i+1} chunks of {chunk_size} rows from the candidate sentences', end='\r')
    print()

    # Concatenate
    TRI_df = pd.concat(TRI_df_list, axis=0)
    nonTRI_df_sample = pd.concat(nonTRI_df_sample_list, axis=0)

    # Save
    TRI_df.to_csv(raw_TRI_path, index=False, sep='\t')
    nonTRI_df_sample.to_csv(raw_nonTRI_sample_path, index=False, sep='\t')

    return

def load_preprocess_df(df_path: str) -> pd.DataFrame:
    df = load_df(df_path)

    df['PMID'] = df['#SentenceID'].apply(lambda row: row.split(':')[1])


    # Change REPRESION to REPRESSION
    df['MoR'] = df['MoR'].replace('REPRESION', 'REPRESSION')

    # Change Valid Score to TRI Score
    df = df.rename(columns={'Valid score': 'TRI score'})

    # 'Mutated_TF' was abandoned as it was in a very small % of sentences and did not give any useful information.
    # df['Mutated_TF'] = np.where(df['Mutated Genes'] != '',
    #                             df.apply(lambda row: row['TF Id'] in row['Mutated Genes'].split(';'), axis=1),
    #                             False)
    # Remove mutated genes column from the df
    df = df.drop(columns=['Mutated Genes', 'Mutation offsets'])

    # Assertions
    assert set(df['MoR'])  == {'ACTIVATION', 'REPRESSION', 'UNDEFINED'}, f'MoR column has unexpected values: {set(df["MoR"])}'
    assert len(df[df['TF'].isnull()]) == 0
    assert len(df[df['TG'].isnull()]) == 0
    assert len(df[df['TF Id'].isnull()]) == 0
    assert len(df[df['TG Id'].isnull()]) == 0

    return df

def remove_duplicates(TRI_df: pd.DataFrame) -> None:
    '''
    When there's more than 1 TF/TG ID combination, in a same sentence, that is labelled as TRI, 
    the combination with the highest prediction from the model will be preserved, while the others will be dropped. 
    '''

    # Ensure colums are updated
    TRI_df['PMID+Sent']         = TRI_df['#SentenceID'].apply(lambda row: row.split(':')[1]+'|'+row.split(':')[4])
    TRI_df['TRI Id']            = TRI_df['TF Id'] + '|' + TRI_df['TG Id']
    TRI_df['PMID+Sent+TRI_Id']  = TRI_df['PMID+Sent'] + '|' + TRI_df['TRI Id']
    TRI_df.drop(columns=['PMID+Sent', 'TRI Id'], inplace=True)

    # Only keep the TRI sentence labelled as TRI with the highest score:
    TRI_df.sort_values(by=['TRI score'], ascending=False, inplace=True)
    TRI_df.drop_duplicates(subset=['PMID+Sent+TRI_Id'], keep="first", inplace=True)
    return

def add_symbols_TaxID(df: pd.DataFrame, EntrezIDtoSymbol_path: str) -> pd.DataFrame:
    ''''''
    with open(EntrezIDtoSymbol_path, 'r') as f:
        EntrezIDtoSymbol = json.load(f)

    def ID_to_(row, var):
        IDs = row.split(";")
        Vars = [EntrezIDtoSymbol.get(ID, {var: ''})[var] for ID in IDs]
        return ";".join(Vars)
    
    def ID_to_name(row):
        return ID_to_(row, 'Name')

    def ID_to_TaxID(row):
        return ID_to_(row, 'TaxID')

    # Drop rows without TG Id 
    m_drop = df[f'TG Id'] == ''
    print(f"{m_drop.sum()} sentences are dropped as their TG is not normalised\n")
    df = df[~m_drop].copy()

    for T in ('TF', 'TG'):
        # Get TF/TG Symbol and TaxID for the rest
        df[f'{T} Symbol'] = df[f'{T} Id'].apply(ID_to_name)
        df[f'{T} TaxID']  = df[f'{T} Id'].apply(ID_to_TaxID)
    
    m = (df['TF Symbol'].isna() | df['TG Symbol'].isna())
    assert m.sum() == 0, "Some gene symbols were not retrieved from EntrezID."

    return df

def remove_other_species(df: pd.DataFrame, TaxID: dict) -> pd.DataFrame:
    '''Remove rows where the TaxID is cross-species'''
    TaxIDs = {'9606', '10090', '10116'}

    for T in ('TF', 'TG'):
        # Keep only 1 copy of each TaxID
        df.loc[:, f'{T} TaxID'] = df.loc[:, f'{T} TaxID'].apply(lambda x: ';'.join(set(x.split(';'))))

        # Only keep rows with 1 TaxID
        df = df[df[f'{T} TaxID'].isin(TaxIDs)]

    return df

# ---------------
# TODO - Remove. This is outdated. Ensure HGNC is not mentioned anywhere
def add_HGNC_symbols(ExTRI2_df: pd.DataFrame, orthologs_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''
    use ortholog dicts in orthologs_path (downloaded from HGNC) to get HGNC orthologs for mouse, rat, & human HGNC IDs
    '''

    ## HELPER FUNCTIONS
    def get_unique_HGNC_symbol_index(row):
        '''
        Helper function to get a unique HGNC symbol for each Entrez ID, when there are multiple options.
        1) Return exact lowercase match, if any
        2) Closest string match (case-insensitive)
        3) If no match, return NaN
        '''

        # If there's only one symbol, return index 0
        if ';' not in row['human_symbol']:
            return 0
        
        # Get all candidates
        candidates = row["human_symbol"].split(";")

        # 1) Return exact lowercase match, if any
        for i, c in enumerate(candidates):
            if c.lower() == row['symbol'].lower():
                return i

        # 2) Closest string match (case-insensitive)
        matches = difflib.get_close_matches(row['symbol'].upper(), [c.upper() for c in candidates], n=1, cutoff=0.0)
        if matches:
            # return the original candidate (not the uppercased one)
            for i, c in enumerate(candidates):
                if c.upper() == matches[0]:
                    return i

        # 3) No matches found - raise error
        raise ValueError(f"No match found for {row['symbol']}")
    
    def assign_unique_fields(row):
        '''Helper function to assign unique fields based on index'''
        idx = get_unique_HGNC_symbol_index(row)

        return pd.Series({
            "unique_human_symbol": row["human_symbol"].split(";")[idx],
            "unique_hgnc_id": row["hgnc_id"].split(";")[idx],
            "unique_human_entrez_gene": row["human_entrez_gene"].split(";")[idx],
        })

    def fill_ortholog_column(id, column):
        '''Helper function to fill ortholog columns'''
        result = []
        for entrez_gene in id.split(";"):
            result.append(orthologs_map[entrez_gene][column]) if entrez_gene in orthologs_map else "-"
        return ";".join(result)

    # Create empty df to store orthologs
    orthologs = pd.DataFrame()

    # Get mouse & rat orthologs
    for rodent in ['mouse', 'rat']:
        # Load as string
        hgnc_df = load_df(orthologs_path + f"human_{rodent}_hcop_fifteen_column.txt")
        hgnc_df = hgnc_df.rename(columns={f"{rodent}_entrez_gene": "entrez_gene", f"{rodent}_symbol": "symbol"})
        hgnc_df['TaxID'] = '10090' if rodent == 'mouse' else '10116'
        orthologs = pd.concat([orthologs, hgnc_df])

    # Get human HNGC symbols
    h_hgnc_df = load_df(orthologs_path + "hgnc_human.tsv")
    h_hgnc_df = h_hgnc_df.rename(columns={"HGNC ID": "hgnc_id", "NCBI Gene ID": "entrez_gene", "Approved symbol": "symbol"})
    h_hgnc_df['human_entrez_gene'] = h_hgnc_df['entrez_gene']
    h_hgnc_df['human_symbol'] = h_hgnc_df['symbol']
    h_hgnc_df['TaxID'] = '9606'
    orthologs = pd.concat([orthologs, h_hgnc_df])

    # Keep only IDs present in the ExTRI2_df
    TF_ids = {j for id in ExTRI2_df['TF Id'].unique() for j in id.split(';')}
    TG_ids = {j for id in ExTRI2_df['TG Id'].unique() for j in id.split(';')}
    TF_TG_ids = TF_ids | TG_ids
    orthologs = orthologs[orthologs['entrez_gene'].isin(TF_TG_ids)]

    # Remove all rows that don't have a symbol, human entrez ID, or hgnc ID
    m = (orthologs['symbol'] != '-') & (orthologs['human_entrez_gene'] != '-') & (orthologs['hgnc_id'] != '-')
    orthologs = orthologs[m]

    # === Resolve cases where entrez_gene has a 1-to-many mapping with symbol. Can be fixed by eliminating all "LOCXXXX" & "GmXXXX" ones. ===
    # Get entrez IDs with multiple symbols
    entrezIDs_with_multiple_symbols = (
        orthologs[['entrez_gene', 'symbol']]
        .drop_duplicates()
        .groupby('entrez_gene')
        .filter(lambda g: len(g) > 1)['entrez_gene']
        .unique()
    )
    # Discard all rows with 'LOCXXXX' or 'GmXXXX' symbols for these entrez IDs (& assert this solves all cases)
    m_to_discard = orthologs['entrez_gene'].isin(entrezIDs_with_multiple_symbols) & (orthologs['symbol'].str.contains(r'^(?:LOC|Gm\d+)', regex=True))
    print(f"Discarding {m_to_discard.sum()} rows with 'LOCXXXX' or 'GmXXXX' symbols: {', '.join(orthologs[m_to_discard]['symbol'].unique())}\n")
    orthologs = orthologs[~m_to_discard]
    assert (orthologs[['entrez_gene', 'symbol']].drop_duplicates()['entrez_gene'].duplicated().sum() == 0), "There are still entrez IDs with multiple symbols after discarding 'LOCXXXX' and 'GmXXXX' ones"


    # Join with ';' when an EntrezID has more than 1 human ortholog
    agg_funcs = {
        "symbol": lambda x: ';'.join(x.unique()),
        "TaxID": lambda x: ';'.join(x.unique()),
        "human_entrez_gene": lambda x: ';'.join(x),
        "hgnc_id": lambda x: ';'.join(x),
        "human_symbol": lambda x: ';'.join(x),
    }
    orthologs = orthologs.groupby(['entrez_gene']).agg(agg_funcs).reset_index()

    # Add NFKB & AP1 orthologs
    for dimer in ['NFKB', 'AP1']:
        orthologs = pd.concat([orthologs, pd.DataFrame([{
            'entrez_gene': f'Complex:{dimer}',
            'symbol': dimer,
            'TaxID': '9606',
            'human_entrez_gene': f'Complex:{dimer}',
            'hgnc_id': f'Complex:{dimer}',
            'human_symbol': dimer,
        }])], ignore_index=True)

    # Add columns with 1-to-1 mapping (unique_human_symbol, unique_hgnc_id, unique_human_entrez_gene), by using either exact match or closest match of symbol wrt human_symbol
    orthologs = orthologs.join(orthologs.apply(assign_unique_fields, axis=1))

    # Show how many orthologs we get
    print(f"We get ortholog info for {len(orthologs)}/{len(TF_TG_ids)} Gene IDs\n")

    # Fill in ortholog columns for TFs and TGs using the "unique" fields for a 1-to-1 mapping
    orthologs_map = orthologs.set_index('entrez_gene').to_dict(orient='index')
    for T in ('TF', 'TG'):
        ExTRI2_df[f"{T}_human_entrez_gene"] = ExTRI2_df[f'{T} Id'].apply(lambda id: fill_ortholog_column(id, "unique_human_entrez_gene"))
        ExTRI2_df[f"{T}_hgnc_id"]           = ExTRI2_df[f'{T} Id'].apply(lambda id: fill_ortholog_column(id, "unique_hgnc_id"))
        ExTRI2_df[f"{T}_human_symbol"]      = ExTRI2_df[f'{T} Id'].apply(lambda id: fill_ortholog_column(id, "unique_human_symbol"))

    return ExTRI2_df, orthologs

def add_human_orthologs(ExTRI2_df: pd.DataFrame, ensembl_folder: str, EntrezID_to_Symbol_p: str) -> pd.DataFrame:
    '''
    Add human orthologs, using Ensembl Genes 115 release.
    Returns: 
        ExTRI2_df with four new columns: TF_human_Id, TG_human_Id, TF_HGNC_ID, TG_HGNC_ID
    '''

    # Load Ensembl orthologs (Ensembl Genes 115 version)
    ensembl_orthologs = pd.read_csv(ensembl_folder + 'mart_export_human_orthologs.txt', sep='\t', dtype='str')
    ensembl_orthologs.rename(columns={'Norway rat - BN/NHsdMcwi gene stable ID': 'Norway rat gene stable ID'}, inplace=True)

    # Only keep relevant rows and columns
    m = ensembl_orthologs['Norway rat gene stable ID'].isna() & ensembl_orthologs['Mouse gene stable ID'].isna()
    ensembl_orthologs = ensembl_orthologs[~m][['Gene stable ID', 'Norway rat gene stable ID', 'Mouse gene stable ID']]

    # Load the mappings from Ensembl Gene IDs to NCBI Gene IDs
    ensembl_human_df = pd.read_csv(ensembl_folder + 'mart_export_homo_sapiens_Ensembl_to_NCBI.txt', sep='\t', dtype='str')
    ensembl_mouse_df = pd.read_csv(ensembl_folder + 'mart_export_mus_musculus_Ensembl_to_NCBI.txt', sep='\t', dtype='str')
    ensembl_rat_df = pd.read_csv(ensembl_folder + 'mart_export_rattus_norvegicus_Ensembl_to_NCBI.txt', sep='\t', dtype='str')

    # Build mapping dictionaries from Ensembl ID to NCBI Gene ID
    ensembl_human_map = dict(zip(ensembl_human_df['Gene stable ID'], ensembl_human_df['NCBI gene (formerly Entrezgene) ID']))
    ensembl_mouse_map = dict(zip(ensembl_mouse_df['Gene stable ID'], ensembl_mouse_df['NCBI gene (formerly Entrezgene) ID']))
    ensembl_rat_map   = dict(zip(ensembl_rat_df['Gene stable ID'],   ensembl_rat_df['NCBI gene (formerly Entrezgene) ID']))

    # Add the NCBI Gene ID to the ortholog dataframe
    ensembl_orthologs['Human NCBI Gene ID'] = ensembl_orthologs['Gene stable ID'].map(ensembl_human_map)
    ensembl_orthologs['Mouse NCBI Gene ID'] = ensembl_orthologs['Mouse gene stable ID'].map(ensembl_mouse_map)
    ensembl_orthologs['Rat NCBI Gene ID']   = ensembl_orthologs['Norway rat gene stable ID'].map(ensembl_rat_map)

    # Build mapping from NCBI Gene ID to Human NCBI Gene ID
    mouse_map = dict(zip(ensembl_orthologs['Mouse NCBI Gene ID'], ensembl_orthologs['Human NCBI Gene ID']))
    rat_map   = dict(zip(ensembl_orthologs['Rat NCBI Gene ID'],   ensembl_orthologs['Human NCBI Gene ID']))
    orthologs_map = {k: v for k, v in (mouse_map | rat_map).items() if pd.notna(k) and pd.notna(v)} | {'Complex:NFKB': 'Complex:NFKB', 'Complex:AP1': 'Complex:AP1'}

    # Build mapping from NCBI Gene ID to HGNC ID
    hgnc_map = dict(zip(ensembl_human_df['NCBI gene (formerly Entrezgene) ID'], ensembl_human_df['HGNC ID']))
    hgnc_map = {k: v for k, v in hgnc_map.items() if pd.notna(k) and pd.notna(v)} | {'Complex:NFKB': 'Complex:NFKB', 'Complex:AP1': 'Complex:AP1'}

    def map_NCBI_ID_to_human_NCBI_ID(orthologs_map, gene_ids, TaxID):
        '''Map gene IDs (human, mouse, rat) to human NCBI Gene IDs, using Ensembl 115 release orthologs.'''
        if TaxID == '9606':
            return gene_ids

        return ";".join([orthologs_map.get(id, '-') for id in gene_ids.split(';')])
        

    def map_NCBI_ID_to_HGNC_ID(hgnc_map, gene_ids):
        '''Map human gene ID to HGNC ID.'''
        human_ids = []
        for id in gene_ids.split(';'):
            human_ids.append(hgnc_map.get(id, '-'))
        return ";".join(human_ids)


    with open(EntrezID_to_Symbol_p, 'r') as f:
        EntrezIDtoSymbol = json.load(f)

    def map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, gene_ids):
        '''Map human gene ID to gene symbol.'''
        human_symbols = [EntrezIDtoSymbol.get(id, {'Name': ''})['Name'] for id in gene_ids.split(';')]
        return ";".join(human_symbols)
        print(gene_ids, human_symbols)

    for T in ('TF', 'TG'):
        # Add Human NCBI ID into the final dataframe
        ExTRI2_df[f'{T}_human_Id'] = ExTRI2_df.apply(lambda row: map_NCBI_ID_to_human_NCBI_ID(orthologs_map, row[f'{T} Id'], row[f'{T} TaxID']), axis=1)

        # Add HGNC ID into the final dataframe
        ExTRI2_df[f'{T}_HGNC_ID'] = ExTRI2_df.apply(lambda row: map_NCBI_ID_to_HGNC_ID(hgnc_map, row[f'{T}_human_Id']), axis=1)
    
        # Add human symbols
        ExTRI2_df[f'{T}_human_symbol'] = ExTRI2_df.apply(lambda row: map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, row[f'{T}_human_Id']), axis=1)
    
    # LOG
    # Show how many orthologs / HGNCs are missing
    m = (ExTRI2_df['TF_human_Id'] == '-') | (ExTRI2_df['TG_human_Id'] == '-')
    print(f"No orthologs found for           {m.sum() / len(ExTRI2_df):.2%} of the rows ({m.sum()} rows)")

    m = ExTRI2_df['TF_human_Id'].str.contains('-') | ExTRI2_df['TG_human_Id'].str.contains('-')
    print(f"Some ortholog is not present in  {m.sum() / len(ExTRI2_df):.2%} of the rows ({m.sum()} rows)")

    m = (ExTRI2_df['TF_HGNC_ID'] == '-') | (ExTRI2_df['TG_HGNC_ID'] == '-')
    print(f"No HGNCs for                     {m.sum() / len(ExTRI2_df):.2%} of the rows ({m.sum()} rows)")

    m = ExTRI2_df['TF_HGNC_ID'].str.contains('-') | ExTRI2_df['TG_HGNC_ID'].str.contains('-')
    print(f"Some HGNC is not present in      {m.sum() / len(ExTRI2_df):.2%} of the rows ({m.sum()} rows)")

    return ExTRI2_df
# ----------------

# --- Orthologs functions

def fetch_entrez_symbols(gene_ids, email="you@example.com", batch=200) -> dict:
    """
    Query NCBI's Entrez Gene (E-utilities) API to get gene symbols for a list of Entrez IDs.

    Returns:
        dict: {gene_id (str): gene_symbol (str)} for successfully retrieved genes.
    """
    geneID_to_symbol = {}
    gene_ids = list(map(str, gene_ids))  # ensure strings for URL joining

    for i in range(0, len(gene_ids), batch):
        # Prepare a batch of IDs
        ids = ",".join(gene_ids[i:i + batch])

        # Call NCBI's esummary endpoint
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "gene", "id": ids, "retmode": "json", "email": email}
        )
        r.raise_for_status()  # raise an exception if HTTP request failed

        # Parse JSON response
        docs = r.json().get("result", {})

        # Extract gene symbols
        for k, v in docs.items():
            if k == "uids":  # skip the metadata field
                continue
            symbol = v.get("name")
            if symbol:       # only add if symbol exists
                geneID_to_symbol[k] = symbol

        # Respect NCBI rate limit (~3 requests per second)
        time.sleep(0.34)

    return geneID_to_symbol

def create_rodent_NCBI_to_human_NCBI_mapping(ensembl_folder: str) -> tuple[dict, dict]:
    '''
    Create a mapping from rodent NCBI IDs to human NCBI IDs, using the Ensembl 115 release data.

    The mapping will be created following this chain:
    rodent NCBI ID > rodent Ensembl ID > human Ensembl ID > human NCBI ID

    Return:
    - rodentID_to_humanID: dict, mapping from rodent NCBI Gene ID to a list of human NCBI Gene IDs
    '''

    # Load the mappings from Ensembl Gene IDs to NCBI Gene IDs
    ensembl_human_df = pd.read_csv(ensembl_folder + 'mart_export_homo_sapiens_Ensembl_to_NCBI.txt', sep='\t', dtype='str')
    ensembl_mouse_df = pd.read_csv(ensembl_folder + 'mart_export_mus_musculus_Ensembl_to_NCBI.txt', sep='\t', dtype='str')
    ensembl_rat_df = pd.read_csv(ensembl_folder + 'mart_export_rattus_norvegicus_Ensembl_to_NCBI.txt', sep='\t', dtype='str')

    ensembl_orthologs_df = pd.read_csv(ensembl_folder + 'mart_export_human_orthologs.txt', sep='\t', dtype='str')
    ensembl_orthologs_df.rename(columns={'Norway rat - BN/NHsdMcwi gene stable ID': 'Norway rat gene stable ID'}, inplace=True)

    # Create rodent ID to TaxID mapping
    rodentID_to_taxID = {}
    for NCBI_ID in ensembl_mouse_df['NCBI gene (formerly Entrezgene) ID'].unique():
        rodentID_to_taxID[NCBI_ID] = 10090
    for NCBI_ID in ensembl_rat_df['NCBI gene (formerly Entrezgene) ID'].unique():
        rodentID_to_taxID[NCBI_ID] = 10116

    # 1) rodent NCBI ID > rodent Ensembl ID
    ensembl_mouse_df = ensembl_mouse_df[ensembl_mouse_df['NCBI gene (formerly Entrezgene) ID'].notna()]
    ensembl_rat_df = ensembl_rat_df[ensembl_rat_df['NCBI gene (formerly Entrezgene) ID'].notna()]
    ensembl_rodent_df = pd.concat([ensembl_mouse_df, ensembl_rat_df], ignore_index=True).drop_duplicates()
    rodent_grouped = (ensembl_rodent_df.groupby('NCBI gene (formerly Entrezgene) ID')['Gene stable ID'].agg(set).reset_index())
    rodent_geneID_to_ensembl = dict(zip(rodent_grouped['NCBI gene (formerly Entrezgene) ID'], rodent_grouped['Gene stable ID']))

    # 2) rodent Ensembl ID > human Ensembl ID
    mouse_grouped = (ensembl_orthologs_df.groupby('Mouse gene stable ID')['Gene stable ID'].agg(set).reset_index())
    rat_grouped = (ensembl_orthologs_df.groupby('Norway rat gene stable ID')['Gene stable ID'].agg(set).reset_index())

    rodent_ensembl_to_human_ensembl = {
        **dict(zip(mouse_grouped['Mouse gene stable ID'], mouse_grouped['Gene stable ID'])), 
        **dict(zip(rat_grouped['Norway rat gene stable ID'], rat_grouped['Gene stable ID']))
    }

    # 3) human Ensembl ID > human NCBI ID
    ensembl_human_df = ensembl_human_df[ensembl_human_df['NCBI gene (formerly Entrezgene) ID'].notna()].drop_duplicates()
    grouped = (ensembl_human_df.groupby('Gene stable ID')['NCBI gene (formerly Entrezgene) ID'].agg(set).reset_index())
    human_ensembl_to_geneID = dict(zip(grouped['Gene stable ID'], grouped['NCBI gene (formerly Entrezgene) ID']))

    # 4) rodent NCBI ID > human NCBI ID (where multiple are joined by ;)
    rodentID_to_humanID = {}
    for rodent_id, rodent_ensembl_ids in rodent_geneID_to_ensembl.items():
        human_ids = set()
        for rodent_ensembl_id in rodent_ensembl_ids:
            human_ensembl_ids = rodent_ensembl_to_human_ensembl.get(rodent_ensembl_id, [''])
            for ensembl_id in human_ensembl_ids:
                human_ids.update(human_ensembl_to_geneID.get(ensembl_id, ['']))
        human_ids.discard('')  # Remove empty strings if any
        rodentID_to_humanID[rodent_id] = ";".join(human_ids) if human_ids else 'None'

    return rodentID_to_humanID, rodentID_to_taxID

def create_human_NCBI_to_HGNC_mapping(ensembl_folder: str) -> dict:
    '''
    Create a mapping from human NCBI Gene IDs to HGNC IDs, using the Ensembl 115 release data.
    If multiple HGNC IDs are found (only in 1% of cases), the NCBI ID is discarded.

    Return:
    - humanID_to_HGNC: dict, mapping from human NCBI Gene ID to a list of HGNC IDs
    '''

    # Load the mappings from Ensembl Gene IDs to NCBI Gene IDs
    ensembl_human_df = pd.read_csv(ensembl_folder + 'mart_export_homo_sapiens_Ensembl_to_NCBI.txt', sep='\t', dtype='str')

    # Filter rows with both HGNC and NCBI IDs, and remove duplicates
    m = ensembl_human_df['HGNC ID'].notna() & ensembl_human_df['NCBI gene (formerly Entrezgene) ID'].notna()
    ensembl_human_df = ensembl_human_df[m][['HGNC ID', 'NCBI gene (formerly Entrezgene) ID']].drop_duplicates()

    # Drop rows with NCBI IDs that have multiple HGNC IDs
    m = ensembl_human_df['NCBI gene (formerly Entrezgene) ID'].duplicated(keep=False)
    print(f"Discarding from the NCBI_Gene_ID-to_HGNC_ID mapping {ensembl_human_df['NCBI gene (formerly Entrezgene) ID'].duplicated().sum()} Human Gene IDs out of {ensembl_human_df['NCBI gene (formerly Entrezgene) ID'].nunique()} as they map to multiple HGNCs")
    ensembl_human_df = ensembl_human_df[~m]

    # Create a one-to-one dictionary
    humanID_to_HGNC = dict(zip(ensembl_human_df['NCBI gene (formerly Entrezgene) ID'], ensembl_human_df['HGNC ID']))

    return humanID_to_HGNC

def create_orthologs_df(ensembl_folder: str, EntrezID_to_Symbol_TRI_p: str, ExTRI2_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Create a dataframe with the orthologs mapping from rodent NCBI Gene ID to human NCBI Gene ID, HGNC ID, and human gene symbol.

    The mapping will be created following this chain:
    rodent NCBI ID > rodent Ensembl ID > human Ensembl ID > human NCBI ID

    Return:
    - orthologs_df: pd.DataFrame, with columns: Gene_ID, human_gene_ID, HGNC_ID, Gene Symbol, human_gene_symbol
     - As well as unique mappings using same name, or first family name rule
    '''

    def map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, gene_ids):
        '''Map human gene ID to gene symbol.'''
        gene_symbols = [EntrezIDtoSymbol.get(id, {'Name': 'None'})['Name'] for id in gene_ids.split(';')]
        return ";".join(gene_symbols)

    def get_unique_human_symbol_index(row):
        '''
        Get the index of the unique human gene symbol following these rules:
        1) Return exact uppercase match, if any
        2) First gene family member (e.g. ACSM2A for ACSM2; if multiple, take the one with the smallest numeric suffix)
        3) If no match, return None
        '''
        # Get rodent symbol & human symbols
        rodent_symbol = row['gene_symbol'].upper()
        human_symbols = row['human_gene_symbol'].upper().split(";")

        # If there's only one symbol, return index 0
        if len(human_symbols) == 1:
            return 0

        # 1) Exact match (case-insensitive)
        for i, c in enumerate(human_symbols):
            if c == rodent_symbol:
                return i

        # 2) Apply first gene family member rule:
        # Split the symbol into stem, number, letter, number2
        m = re.match(r'^([A-Z]+)(\d+)?([A-Z]?)(\d+)?$', rodent_symbol)
        stem = m.group(1) if m else rodent_symbol
        num = m.group(2) if (m and m.group(2)) else None
        letter = m.group(3) if (m and m.group(3)) else None
        num2 = m.group(4) if (m and m.group(4)) else None

        # Get human symbols that start with the same stem
        candidates = [(i, hs) for i, hs in enumerate(human_symbols) if hs.startswith(stem)]
        if not candidates:
            return None
        
        # Prefer candidates that match the stem
        exact = [i for i, hs in candidates if hs == stem]
        if exact:
            return exact[0]
            
        # If there is a number, prefer same numeric family (e.g. ADH1A over ADH2A)
        if num:
            same_family = [(i, hs) for i, hs in candidates if hs.startswith(stem + num)]

            if same_family:
                # Prefer exact stem+num (if present), else smallest suffix
                exact = [i for i, hs in same_family if hs == stem + num]
                if exact:
                    return exact[0]
                
                same_family.sort(key=lambda x: x[1])
                return same_family[0][0]
            
        # Otherwise, just pick smallest lexicographic suffix (first family member)
        candidates.sort(key=lambda x: x[1])
        return candidates[0][0]

    def assign_unique_human_fields(row):
        '''Helper function to assign unique fields based on index'''
        idx = get_unique_human_symbol_index(row)

        if idx is None:
            return pd.Series({
                'unique_human_gene_ID': 'None',
                'unique_human_gene_symbol':  'None',
            })

        return pd.Series({
            'unique_human_gene_ID': row['human_gene_ID'].split(';')[idx],
            'unique_human_gene_symbol': row['human_gene_symbol'].split(';')[idx],
        })

    # Create mappings
    rodentID_to_humanIDs, rodentID_to_taxID = create_rodent_NCBI_to_human_NCBI_mapping(ensembl_folder)
    humanID_to_HGNC = create_human_NCBI_to_HGNC_mapping(ensembl_folder)

    # Load EntrezID to Symbol mapping
    with open(EntrezID_to_Symbol_TRI_p, 'r') as f:
        EntrezIDtoSymbol = json.load(f)

    # Create a dataframe with the orthologs mapping
    orthologs_df = pd.DataFrame.from_dict(rodentID_to_humanIDs, orient='index').reset_index()
    orthologs_df.columns = ['Gene_ID', 'human_gene_ID']

    # Only keep rows that are present in the ExTRI2_df
    geneIDs_in_ExTRI2 = {id for col in ['TF Id', 'TG Id'] for ids in ExTRI2_df[col].unique() for id in ids.split(';')}
    orthologs_df = orthologs_df[orthologs_df['Gene_ID'].isin(geneIDs_in_ExTRI2)].reset_index(drop=True)

    # EntrezIDtoSymbol only contains IDs in ExTRI. Retrieve missing ones from Entrez
    all_humanIDs = {id  for ids in orthologs_df['human_gene_ID'].unique() for id in ids.split(';')}
    geneIDs_w_missing_symbol = [id for id in all_humanIDs if (id not in EntrezIDtoSymbol) and (id not in ['-', ''])]
    print(f"Fetching {len(geneIDs_w_missing_symbol)} missing human gene symbols from Entrez...")
    if geneIDs_w_missing_symbol:
        # Fetch missing symbols from Entrez
        missing_symbols = fetch_entrez_symbols(geneIDs_w_missing_symbol)
        EntrezIDtoSymbol.update({id: {'Name': symbol, 'TaxID': '9606'} for id, symbol in missing_symbols.items()})

    # Add TaxID
    orthologs_df['TaxID'] = orthologs_df['Gene_ID'].apply(lambda id: rodentID_to_taxID[id])

    # Add both rodent and human gene symbol
    orthologs_df['gene_symbol'] = orthologs_df['Gene_ID'].apply(lambda id: map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, id))
    orthologs_df['human_gene_symbol'] = orthologs_df['human_gene_ID'].apply(lambda id: map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, id))

    # Add unique human gene ID and symbol
    orthologs_df = orthologs_df.join(orthologs_df.apply(assign_unique_human_fields, axis=1))

    # Add the human IDs to the orthologs df for completeness
    human_df = pd.DataFrame(list(humanID_to_HGNC.items()), columns=['Gene_ID', 'HGNC_ID'])
    human_df = human_df[human_df['Gene_ID'].isin(geneIDs_in_ExTRI2)].reset_index(drop=True) # Drop rows not in ExTRI2
    human_df['gene_symbol'] = human_df['Gene_ID'].apply(lambda id: map_NCBI_ID_to_Symbol(EntrezIDtoSymbol, id))
    human_df['TaxID'] = '9606'

    # For consistency with orthologs_df columns:
    human_df['human_gene_ID'] = human_df['Gene_ID']
    human_df['human_gene_symbol'] = human_df['gene_symbol']
    human_df['unique_human_gene_ID'] = human_df['Gene_ID']
    human_df['unique_human_gene_symbol'] = human_df['gene_symbol']

    # Join with orthologs_df
    orthologs_df = pd.concat([orthologs_df, human_df], ignore_index=True)

    # Add HGNC symbol
    orthologs_df['HGNC_ID'] = orthologs_df['human_gene_ID'].apply(lambda human_gene_ID: ';'.join(humanID_to_HGNC.get(id, 'None') for id in human_gene_ID.split(';')))
    orthologs_df['unique_HGNC_ID'] = orthologs_df['unique_human_gene_ID'].apply(lambda human_gene_ID: humanID_to_HGNC.get(human_gene_ID, 'None'))

    # Add API & NFKB complexes:
    orthologs_df = pd.concat([orthologs_df, pd.DataFrame([
        {'Gene_ID': 'Complex:NFKB', 'human_gene_ID': 'Complex:NFKB', 'gene_symbol': 'Complex:NFKB', 'human_gene_symbol': 'Complex:NFKB', 'unique_human_gene_ID': 'Complex:NFKB', 'unique_human_gene_symbol': 'Complex:NFKB', 'HGNC_ID': 'Complex:NFKB', 'unique_HGNC_ID': 'Complex:NFKB'},
        {'Gene_ID': 'Complex:AP1',  'human_gene_ID': 'Complex:AP1',  'gene_symbol': 'Complex:AP1',  'human_gene_symbol': 'Complex:AP1',  'unique_human_gene_ID': 'Complex:AP1',  'unique_human_gene_symbol': 'Complex:AP1',  'HGNC_ID': 'Complex:AP1',  'unique_HGNC_ID': 'Complex:AP1'},
    ])], ignore_index=True)

    return orthologs_df

def add_ortholog_columns(ExTRI2_df: pd.DataFrame, orthologs_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Apply unique orthologs mapping to the ExTRI2 dataframe.
    '''
    # Get orthologs_map from the orthologs_df
    orthologs_map = {
        row['Gene_ID']: {
            'human_Id': row['unique_human_gene_ID'],
            'human_symbol': row['unique_human_gene_symbol'],
            'HGNC_Id': row['unique_HGNC_ID'],
            # 'multiple_human_Id': row['human_gene_ID'],
            # 'multiple_human_symbol': row['human_gene_symbol'],
            # 'multiple_HGNC_Id': row['HGNC_ID']
        }
        for _, row in orthologs_df.iterrows()
    }

    # Apply the map to ExTRI2
    for T in ['TF', 'TG']:
        for field in ['human_Id', 'human_symbol', 'HGNC_Id', ]: #'multiple_human_Id', 'multiple_human_symbol', 'multiple_HGNC_Id']:
            ExTRI2_df[f'{T}_{field}'] = ExTRI2_df[f'{T} Id'].apply(
                lambda ids: ';'.join(orthologs_map.get(id, {}).get(field, 'None') for id in ids.split(';') if id.strip())
            )

    # LOG - Show how many orthologs / HGNCs are missing
    for text, field in [('orthologs', 'human_Id'), ('HGNCs', 'HGNC_Id'), ('human symbols', 'human_symbol')]:
        m1 = (ExTRI2_df['TF_' + field] == 'None') | (ExTRI2_df['TG_' + field] == 'None')
        m2 = ExTRI2_df['TF_' + field].str.contains('None') | ExTRI2_df['TG_' + field].str.contains('None')
        print(f"No {text+' in':<17} {m1.sum()} rows ({m1.mean():.2%}). Some missing in {m2.sum()} rows ({m2.mean():.2%})")

    return ExTRI2_df
# ----------------

def add_TF_type(ExTRI2_df: pd.DataFrame, config: dict) -> None:
    '''Assign a TF type to each sentence based on the IDs in the TF_path. If not found, write "-"'''

    def load_TF_set(TF_type):
        TF_path = f"{config['TF_ids_p']}{TF_type}_entrez_code.list"
        TF_set = set()
        with open(TF_path, "r") as f:
            for l in f:
                TF_set.add(l.strip("\n"))
        return TF_set
    
    # Get Dictionary of TF ids sets: 
    TF_sets = {}
    for TF_type in config['TF_types']:
        # Get set of TF ids of that type
        TF_set = load_TF_set(TF_type)
        TF_sets[TF_type] = TF_set
    
        # Fill in the TF type when there's only 1 TF Id
        m = ExTRI2_df['TF Id'].isin(TF_set)
        ExTRI2_df.loc[m, 'TF_type'] = TF_type

    # Remove ll_coTF from the coTF set
    TF_sets['coTF'].difference_update(TF_sets['ll_coTF'])

    # Fill in TF type when there's more than 1 TF Id (e.g. 111654;15977)
    def get_TF_type_list(TF_IDs):
        '''Return the TF types of the TF IDs joined by ;'''
        TF_type_list = [next((TF_type for TF_type, TF_set in TF_sets.items() if id in TF_set), '-') for id in TF_IDs.split(";")]
        return ";".join(TF_type_list)
    m = ExTRI2_df['TF_type'].isna()
    ExTRI2_df.loc[m, 'TF_type'] = ExTRI2_df[m]['TF Id'].apply(get_TF_type_list)

    # Assert there are no TFs without type
    assert (isna_sum := ExTRI2_df['TF_type'].isna().sum()) == 0, f"There are {isna_sum} TFs without type"

    return

def drop_GTFs(ExTRI2_df: pd.DataFrame) -> pd.DataFrame:
    '''Drop all sentences where the TF uniquely maps to a GTF'''
    m_drop = ExTRI2_df["TF_type"].apply(lambda x: set(x.split(';')).issubset({'GTF', '-'}))
    return ExTRI2_df[~m_drop]

def save_df(df, output_p):
    df.to_csv(output_p, index=False, sep='\t')

def postprocess(ExTRI2_df: pd.DataFrame, TRI_sents: bool, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Add metadata, filter out sentences and renormalize entities with common mistakes

    Args:
        TRI_sents: bool, indicates whether to consider sentences labelled as TRI or "not TRI"
        config: dictionary with paths
    '''

    df_type = 'TRI' if TRI_sents else 'nonTRI'
    print(f'### POSTPROCESSING {df_type}_df')

    # Retrieve Symbol & TaxID from Entrez
    save_Symbol_TaxID_dict(ExTRI2_df, config[f'EntrezID_to_Symbol_{df_type}_p'])

    # Filter & add metadata
    if TRI_sents:
        remove_duplicates(ExTRI2_df)
    ExTRI2_df = add_symbols_TaxID(ExTRI2_df, config[f'EntrezID_to_Symbol_{df_type}_p'])
    add_TF_type(ExTRI2_df, config)
    ExTRI2_df = drop_GTFs(ExTRI2_df)
    ExTRI2_df = remove_other_species(ExTRI2_df, TaxID)

    # Fix AP1 & NFKB normalisations
    ExTRI2_df = fix_NFKB_AP1(ExTRI2_df, config)

    # Renormalize & discard sentences with common mistakes
    if TRI_sents:
        # Save version before renormalisation
        save_df(ExTRI2_df, config['TRI_pre_renorm_p'])

        # Renormalise & discard
        renormalize(ExTRI2_df, renormalized_sents_path=config['renormalized_sents_p'])
        ExTRI2_df = discard(ExTRI2_df, discarded_sents_path=config['discarded_sents_p'])
    else:
        renormalize(ExTRI2_df, renormalized_sents_path = None)
        ExTRI2_df = discard(ExTRI2_df, discarded_sents_path = None)

    # Remove duplicates (if formed after renormalization)
    if TRI_sents:
        remove_duplicates(ExTRI2_df)

    # TODO - Remove these two, they're outdated
    # Add HGNC symbols
    # ExTRI2_df, orthologs_df = add_HGNC_symbols(ExTRI2_df, config['orthologs_p'])
    # Add human orthologs
    # ExTRI2_df = add_human_orthologs(ExTRI2_df, config['ensembl_folder'], config[f'EntrezID_to_Symbol_{df_type}_p'])

    # Add orthologs using Ensembl 115 release
    orthologs_df = create_orthologs_df(config['ensembl_folder'], config[f'EntrezID_to_Symbol_{df_type}_p'], ExTRI2_df)
    ExTRI2_df = add_ortholog_columns(ExTRI2_df, orthologs_df)
    print()

    return ExTRI2_df, orthologs_df


def main():
    config = load_config()

    # Only run if not run already
    if not os.path.isfile(config['raw_TRI_p']):
        # Get manageable datasets from the raw ExTRI file (which has millions of candidate sentences)
        get_TRI_nonTRI_raw_dbs(config['raw_results_p'], config['raw_TRI_p'], config['raw_nonTRI_sample_p'], config['chunk_size'])

    # Load raw dataframes
    TRI_df            = load_preprocess_df(config['raw_TRI_p'])
    not_TRI_sample_df = load_preprocess_df(config['raw_nonTRI_sample_p'])

    # Postprocess
    TRI_df, orthologs_df = postprocess(TRI_df,            TRI_sents=True,  config=config)
    not_TRI_sample_df, _ = postprocess(not_TRI_sample_df, TRI_sents=False, config=config)

    # Save
    save_df(TRI_df, config['final_ExTRI2_p'])
    save_df(not_TRI_sample_df, config['nonTRI_sample_p'])
    
    # Save orthologs_df for reference
    # TODO - This must be updated in Drive
    orthologs_df.to_csv(config['orthologs_final_p'], sep='\t', index=False)

if __name__ == "__main__":
    main()


