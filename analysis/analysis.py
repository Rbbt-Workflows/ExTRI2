'''Functions related to the analysis of the results'''

import pandas as pd
import numpy as np
from IPython.display import display, HTML
import matplotlib.pyplot as plt





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
    config['validated_sents_p']     = RESULTS_P + 'validated_sentences.tsv'

    # TF lists
    config['dbTFs_p']       = DATA_P + 'postprocessing/dbTF_entrez_code.list'
    config['coTFs_p']       = DATA_P + 'postprocessing/coTF_entrez_code.list'
    config['ll_coTFs_p']    = DATA_P + 'postprocessing/ll_coTF_entrez_code.list'
    config['TFs_p']         = DATA_P + 'tf_entrez_code.list'

    return config

def load_TF_set(path: str) -> set:
    with open(path) as f:
        TF_set = {l.strip('\n') for l in f}
    return TF_set

def main():
    config = load_config()

if __name__ == "__main__":

    main()