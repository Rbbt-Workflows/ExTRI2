# Content
This folder contains all original and intermediate datasets required for creating the final ExTRI2 dataset. It only excludes the datasets and models used for training the classifiers, so that the `classifiers_training/` model can be run independently.

The folders are:
* `dataset_improvement/` contains the reannotated tsvs with sentences from the training dataset, as created by `scripts/classifiers_training/prepare_reannotation_Excels.ipynb`
* `external/` contains all datasets that come from outside ExTRI: ortholog lists, TF IDs, and the original training dataset.
* `postprocessing/` contains datasets useful for postprocessing. Specifically
* `pubtator/` contains the abstracts used to obtain the ExTRI2 pipeline.
* `tf_entrez_code.list` has a list of all IDs considered as TFs by the ExTRI2 pipeline.


終 TODO - Explan how each of the datasets has been obtained

## Obtention of `external/` datasets
### human_HGNC_orthologs
This folder contains 3 tsv files with orthology mappings between mouse, rat and human. 

human_mouse & human_rat mappings were obtained from https://www.genenames.org/tools/hcop/. From there, I downloaded (01/07/2024):
* human - mouse ortholog data as a 15 column tab delimited text file
* human - rat ortholog data as a 15 column tab delimited text file
hgnc_human was downloaded from https://www.genenames.org/download/custom/ (25/08/2024) to get the HGNC IDs from human NCBI IDs

> TODO - Be more detailed

### TF_id
> TODO - Explain how I obtained the QuickGO & TFCheckpoint annotations
>

### tsv files
`original_tri_sentences.tsv` is the original tsv file that was used for training the classifiers (see `classifiers_training/`)

`NTNU_extended.tsv` is where `original_tri_sentences.tsv` comes from, after applying a group of filters

> TODO - Should we explain the filters used to obtain the original_train_sentences? We didn't use them all and should explain why right?