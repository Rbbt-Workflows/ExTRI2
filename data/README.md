# Overview
This folder contains lightweight reference files required for the ExTRI2 pipeline. Large and intermediate datasets are not included in the repository and must be generated or downloaded separately. 


This folder contains all original and intermediate datasets required for creating the final ExTRI2 dataset, except for those used for classifiers training (found in `classifiers_training/`), so that classifier training can be done independently from the rest.

# Structure
- `tf_entrez_code.list` contains the list of all IDs considered as TFs by the ExTRI2 pipeline.
- `preprocessing/` contains small preprocessing reference files (e.g. `pubtator_gene_ids.txt`) used to prepare datasets for analysis.
- `postprocessing/` contains configuration or lightweight files used by the postprocessing scripts.
- `external/` contains small reference datasets required for the pipeline (see below for details).

# External datasets
The `external/` folder includes reference datasets required to run the pipeline:
- `TF_id/`: TF IDs from QuickGO and TFCheckpoint  (generated via `scripts/preprocessing/get_NCBI_TF_IDs.ipynb`)
- `ensembl_Release_115_orthologs/`: orthology mappings between mouse, rat, and human (Ensembl release 115, downloaded via BioMart)
- `human_HGNC_orthologs/`: legacy orthology mappings from https://www.genenames.org/tools/hcop/ (kept for reference)
- `all_human_TGs.tsv`: list of all human genes (downloaded from NCBI Datasets, https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9606/, 15/11/2024)

# Data availability
Large datasets (e.g. PubTator abstracts and intermediate pipeline outputs) are not stored in this repository.

To reproduce the results:
1. Download PubTator3 data:
   https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/
2. Run preprocessing scripts:
   - `scripts/preprocessing/get_all_pubtators.sh`
   - `scripts/preprocessing/prepare_pubtator_for_ExTRI2.py`
3. Run the ExTRI2 workflow and postprocessing as described in the main README.
4. Final ExTRI2 dataset and training dataset is available via Zenodo: DOI 10.5281/zenodo.19816074