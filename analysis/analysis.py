'''Functions related to the analysis of the results'''

import pandas as pd
import numpy as np
from IPython.display import display, HTML
import matplotlib.pyplot as plt
from Bio import Entrez

def load_config() -> dict:
    config = {}

    # Folders
    RESULTS_P       = '../results/'
    DATA_P          = '../data/'
    CLASSIFIERS_P   = '../classifiers_training/data/'
    TEMP_P          = DATA_P + 'tmp/'

    # Results
    config['raw_ExTRI2_p']          = RESULTS_P + 'ExTRI2.tsv'
    config['final_ExTRI2_p']        = RESULTS_P + 'ExTRI2_final_resource.tsv'
    config['final_validated_p']     = RESULTS_P + 'validated_sentences.tsv'

    # Classifiers
    config['models_p'] = '../classifiers_training/saved_models/'

    # Other tables
    config['discarded_sents_p']    = DATA_P + 'postprocessing/tables/discarded_sents.tsv'
    config['nfkb_ap1_discarded_p'] = DATA_P + 'postprocessing/tables/NFKB_AP1_discarded_sents.tsv'
    config['all_human_TGs_p']      = DATA_P + 'external/all_human_genes.tsv'
    config['all_TFs_p']            = 'tables/all_TFs.tsv'


    # Tables included in the paper
    config['paper_tables_p'] = DATA_P + 'paper_tables/'
    config['all_considered_TFs_p']  = config['paper_tables_p'] + 'all_considered_TFs.tsv'
    config['TFs_in_ExTRI2_p']       = config['paper_tables_p'] + 'TFs_in_ExTRI2.tsv'
    config['paper_ExTRI2_p']        = config['paper_tables_p'] + 'ExTRI2_final_resource.tsv.gz'
    config['paper_validated_p']     = config['paper_tables_p'] + 'validated_sentences.tsv'
    config['NTNU_dataset_p']        = config['paper_tables_p'] + 'NTNU_training_dataset.tsv'
    config['all_discarded_sents_p'] = config['paper_tables_p'] + 'discarded_sents.tsv'

    return config

def load_TF_set(path: str) -> set:
    with open(path) as f:
        TF_set = {l.strip('\n') for l in f}
    return TF_set

def retrieve_annotations_entrez(id_list):
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)
    webEnv = result["WebEnv"] # type: ignore
    queryKey = result["QueryKey"] # type: ignore
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
    annotationsSummary = annotations['DocumentSummarySet']['DocumentSummary'] # type: ignore

    assert len(id_list) == len(annotationsSummary), f"id_list and annotationsSummary are of different length: {len(id_list)} != {len(annotationsSummary)}"

    return annotationsSummary

def modify_TF_type(df, TFs_df):
    '''
    Modify the "TF_type" column of df by the new categorisation present in TFs_df.
    Rows with no TF type are marked as "-" and then removed.
    '''

    # Build the mapping dict
    id_to_type = dict(zip(
        TFs_df['Gene ID'].astype(str),
        TFs_df['updated TF type'].fillna('-').astype(str)
    ))
    # Add NFKB and AP1 complexes
    id_to_type['Complex:NFKB'] = 'dbTF'
    id_to_type['Complex:AP1'] = 'dbTF'

    # Make sure all unclassified genes are marked with "-"
    assert (set(id_to_type.values()) == {'dbTF', 'coTF', 'coTF candidate', '-'}), "Unexpected TF types found in the mapping dict"

    # Function to map IDs
    def map_ids(row):
        return ";".join(
            id_to_type.get(x.strip(), '-')   # use "-" if unknown
            for x in row['TF Id'].split(";")
        )

    df['new_TF_type'] = df.apply(map_ids, axis=1)

    # Show the difference between old and new TF types
    print("Difference between old and new TF types:")
    display(df[['TF_type', 'new_TF_type']].value_counts(dropna=False).head(5))

    # Replace old TF_type by new_TF_type
    df.drop(columns=['TF_type'], inplace=True)
    df.rename(columns={'new_TF_type': 'TF_type'}, inplace=True)
    # Remove rows with no TF type (all "-")
    n_before = len(df)
    df = df[df['TF_type'].str.split(";").apply(lambda x: set(x) != {'-'})].copy()
    n_after = len(df)
    print(f"Removed {n_before - n_after} ({(n_before - n_after) / n_before:.2%}) rows with no TF type ({n_before} -> {n_after})")

    return df

def save_paper_validated_df(validated_df: pd.DataFrame, TFs_df: pd.DataFrame, path: str) -> pd.DataFrame:
    '''Get the validated_df in the format specified in the paper'''

    # --- Modify TF type to the updated version
    # Drop validated_df repeated column
    assert (validated_df['TF_type_validated'] == validated_df['TF_type']).all(), "The TF types in the validated df have been modified!"
    validated_df.drop(columns=['TF_type_validated'], inplace=True)

    # Modify the TF_type to the updated version
    validated_df = modify_TF_type(validated_df, TFs_df)

    # Only keep relevant columns
    columns = ['#SentenceID', 'Sentence', 'Label', 'MoR', 'Valid?', 'true_label', 'true_MoR', 'Valid score', 'MoR scores',  # TRI - MoR validation
               'TF', 'TF Id', 'TF Symbol', 'TF TaxID',  'TF_is_incorrect', 'TF_correct_mention', 'TF offset',               # TF validation
               'TG',  'TG Id', 'TG Symbol', 'TG TaxID', 'TG_is_incorrect', 'TG_correct_mention', 'Gene offset',             # TG validaiton
               'method', 'renormalisation', 'pre-post', 'TF_type',                                                 # Metadata
               'Other issues', 'Explanation'
               ]
    
    # Only keep relevant columns
    paper_val_df = validated_df[columns].copy()

    # Change column names to match with paper
    paper_val_df.rename(columns={'Label': 'TRI', 'true_label': 'true_TRI', 'Valid score': 'TRI score'}, inplace=True)

    # Change TRI labels to match with paper
    paper_val_df['TRI'] = paper_val_df['TRI'].replace({'FALSE': 'Not TRI', 'TRUE': 'TRI'})
    paper_val_df['true_TRI'] = paper_val_df['true_TRI'].replace({'FALSE': 'Not TRI', 'TRUE': 'TRI'})    
    
    # Save the paper version
    pd.DataFrame.to_csv(paper_val_df, path, sep='\t', index=False, header=True)
    
    return paper_val_df

def combine_discarded_sents(config: dict) -> None:
    '''Combine the 2 tables regarding discarded sentences'''

    discarded_sents_df    = pd.read_csv(config['discarded_sents_p'], sep="\t", dtype=str)
    nfkb_ap1_discarded_df = pd.read_csv(config['nfkb_ap1_discarded_p'], sep="\t", dtype=str)

    discarded_sents_df = pd.concat([discarded_sents_df, nfkb_ap1_discarded_df], ignore_index=True)
    discarded_sents_df.drop(columns=['Unnamed: 0', 'PMID+Sent+TRI_Id', 'TF_type'], inplace=True)  # drop the index column if present

    print("Saved all discarded sentences in", config['all_discarded_sents_p'])
    discarded_sents_df.to_csv(config['all_discarded_sents_p'], sep="\t", index=False)

    return

def create_paper_TF_tables(ExTRI2_df: pd.DataFrame, TFs_df: pd.DataFrame, config: dict) -> pd.DataFrame:


    # Create 'human_TF_type' column, which prioritises the human TF type, and otherwise, it follows the priority order:
    TFtype_priority = {'dbTF': 3, 'coTF': 2, 'coTF candidate': 1} # dbTF > coTF > coTF candidate

    # Function to resolve TF type per human_gene_ID
    def resolve_human_tf_type(df):
        # If only one unique TF type, keep it
        if df['TF type'].nunique() == 1:
            return df['TF type'].iloc[0]
        
        # Prefer the human ortholog (TaxID == 9606)
        human_rows = df[df['TaxID'].astype(str).str.contains('9606')]
        if not human_rows.empty:
            return human_rows['TF type'].iloc[0]
        
        # Otherwise, choose by priority
        tf_types = sorted(df['TF type'].unique(), key=lambda t: TFtype_priority.get(t, 0), reverse=True)
        return tf_types[0]

    # Get TFs that appear in ExTRI2 at least once
    geneIDs_in_ExTRI2 = {id for col in ['TF Id', 'TG Id'] for ids in ExTRI2_df[col].unique() for id in ids.split(';')}
    TFs_df["In ExTRI2"] = TFs_df['Gene ID'].isin(geneIDs_in_ExTRI2)

    # Save all considered TFs (those with a defined TF type)
    m = TFs_df['TF type'].isin(TFtype_priority.keys())
    print(f'We will discard {(~m).sum()} TFs that are not classified into any TF type')
    TFs_df[~m][['Gene ID', 'Symbol', 'TaxID']].to_csv(config['all_considered_TFs_p'], sep="\t", index=False)

    # Create table with only TFs in ExTRI2
    TFs_in_ExTRI2 = (TFs_df
        .drop(columns=['TF type']).rename(columns={'updated TF type': 'TF type'})
        .loc[TFs_df['In ExTRI2'] & (TFs_df['TF type'].isin(TFtype_priority.keys()))].drop(columns=['In ExTRI2'])
    )

    # Create mapping from human_gene_ID → resolved TF type and add it to the main DataFrame
    tf_type_map = (
        TFs_in_ExTRI2.groupby(['human_gene_ID'])[['TF type', 'TaxID']]
        .apply(resolve_human_tf_type).to_dict()
    )
    TFs_in_ExTRI2['human_TF_type'] = TFs_in_ExTRI2['human_gene_ID'].map(tf_type_map)
    TFs_df['human_TF_type'] = TFs_df['human_gene_ID'].map(tf_type_map)

    # Save the table
    TFs_in_ExTRI2.to_csv(config['TFs_in_ExTRI2_p'], sep="\t", index=False)

    return TFs_df

def main():
    config = load_config()

if __name__ == "__main__":

    main()