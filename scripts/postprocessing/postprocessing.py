'''Run from postprocessing_checkings.ipynb'''

import pandas as pd
import numpy as np
import os
import json
from Bio import Entrez
import time
import difflib

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
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
    annotationsSummary = annotations['DocumentSummarySet']['DocumentSummary']

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
    m = df['TF TaxID'].apply(lambda x: ";".join(set(x.split(';')))).isin(TaxIDs)
    m &= df['TG TaxID'].apply(lambda x: ";".join(set(x.split(';')))).isin(TaxIDs)
    
    return df[m]

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

    # Add HGNC symbols
    ExTRI2_df, orthologs_df = add_HGNC_symbols(ExTRI2_df, config['orthologs_p'])

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
    orthologs_df.to_csv(config['orthologs_final_p'], sep='\t', index=False)

if __name__ == "__main__":
    main()


