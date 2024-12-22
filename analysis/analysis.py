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

    # TF lists
    config['dbTFs_p']       = DATA_P + 'postprocessing/dbTF_entrez_code.list'
    config['coTFs_p']       = DATA_P + 'postprocessing/coTF_entrez_code.list'
    config['ll_coTFs_p']    = DATA_P + 'postprocessing/ll_coTF_entrez_code.list'
    config['TFs_p']         = DATA_P + 'tf_entrez_code.list'
    config['all_human_TGs_p'] = DATA_P + 'external/all_human_genes.tsv'

    return config

def load_TF_set(path: str) -> set:
    with open(path) as f:
        TF_set = {l.strip('\n') for l in f}
    return TF_set

# TODO - Remove this if not needed
# def add_extra_cols(ExTRI2_df, validated_df):
#     '''Add extra columns for analysis purposes'''
#     ExTRI2_df['PMID+Sent']          = ExTRI2_df['#SentenceID'].apply(lambda row: row.split(':')[1]+'|'+row.split(':')[4])
#     ExTRI2_df['TRI Id']             = ExTRI2_df['TF Id'] + '|' + ExTRI2_df['TG Id']
#     ExTRI2_df['PMID+Sent+TRI_Id']   = ExTRI2_df['PMID+Sent'] + '|' + ExTRI2_df['TRI Id']
#     ExTRI2_df['PMID+Sent+TRI']      = ExTRI2_df['PMID+Sent'] + '|' + ExTRI2_df['TF'] + '|' + ExTRI2_df['TG']

#     validated_df['PMID+Sent'] = validated_df['#SentenceID'].apply(lambda row: row.split(':')[0]+'|'+row.split(':')[1])
#     validated_df['PMID+Sent+TRI_Id'] = validated_df['PMID+Sent'] + '|' + validated_df['TF Id'] + '|' + validated_df['TG Id']



def retrieve_annotations_entrez(id_list):
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
    annotationsSummary = annotations['DocumentSummarySet']['DocumentSummary']

    assert len(id_list) == len(annotationsSummary), f"id_list and annotationsSummary are of different length: {len(id_list)} != {len(annotationsSummary)}"

    return annotationsSummary

def main():
    config = load_config()

if __name__ == "__main__":

    main()