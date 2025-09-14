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

    # Gene IDs
    config['all_TFs_p']         = 'tables/all_TFs.tsv'
    config['all_human_TGs_p'] = DATA_P + 'external/all_human_genes.tsv'

    # Classifiers
    config['models_p'] = '../classifiers_training/saved_models/'

    # Tables included in the paper
    config['paper_tables_p'] = DATA_P + 'paper_tables/'
    config['all_considered_TFs_p']  = config['paper_tables_p'] + 'all_considered_TFs.tsv'
    config['TFs_in_ExTRI2_p']       = config['paper_tables_p'] + 'TFs_in_ExTRI2.tsv'
    config['paper_ExTRI2_p']        = config['paper_tables_p'] + 'ExTRI2_final_resource.tsv'
    config['paper_validated_p']     = config['paper_tables_p'] + 'validated_sentences.tsv'

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