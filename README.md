# ExTRI2

## Overview
This repository contains the scripts and datasets used for the development of the ExTRI2 pipeline.

## Setup
- [ ] Explain how to setup RBBT

## Structure
- `workflow.rb`. Main script to obtain TRI sentences from a folder of PubTator files.
- `classifiers_training/`. Standalone folder used to obtain the TRI and MoR classifiers to use in `workflow.rb`. See the `README` inside the folder for a more detailed explanation.
- `scripts/` All other scripts to aid the main `worklow.rb` one, including:
  - **Postprocessing:** 
    - `postprocessing/` to prepare the input for the main script, 
    - `classifiers_training` to prepare the data to train the classifiers.
  - **Preprocessing:**
    - `proprocessing` to convert the output to the final ExTRI2 dataset
    - `validation` to prepare the validation sets to manually validate
- `data/` All raw and intermediate files required to run the workflow. See more information in the `README` inside the folder.
- `results/` contains the raw and final ExTRI2 resource, and the validated sentences.
- `analysis/` contains all analysis of the ExTRI2 dataset used for the ExTRI2 paper


## Steps to obtain the results:
The final ExTRI2 dataset required training the classifiers and improving the training dataset, preparing files for the `workflow.rb`, postprocessing the resulting files, and preparing sentences for validation. This was achieved by running the following scripts:

- **Classifiers training:** 
    - Classifiers were trained inside the folder `classifiers_training/` (see `README` there). 
    - The best performing model was chosen with `analysis/classifiers_comparison.ipynb`.
    - Classifier outputs were used to retroactively detect sentences to revise and improve the training dataset. Inside `scripts/classifiers_training`:
        - `prepare_reannotation_Excels.ipynb` prepares the sentences to revise.
        - `update_tri_sentences.ipynb` and `make_train_data.ipynb` update the datasets and prepare the files used for training the models.
- **Preprocessing:**
  - `preprocessing/get_NCBI_TF_IDs.ipynb` obtains the `tf_entrez_code.list` which determines which Gene IDs are considered as TFs
  - `preprocessing/prepare_pubtator_for_ExTRI2.ipynb`



## Preprocessing
Obtaining the files and models required to run the main ExTRI2 workflow: 
1. `data/tf_entrez_code.list` a list of all dbTFs, coTFs & GTFs. Obtained through running `scripts/preprocessing/get_NCBI_TF_IDs.ipynb`.
2. `data/pubtator/` with all PubMed abstracts containing TFs from the above list, in Pubtator format. To obtain, run:

```
cd scripts/preprocessing/
./get_all_pubtators.sh
python prepare_pubtator_for_ExTRI2.py
```

3. `TRI_classifier` and `MoR_classifier` models to classify the sentences. 
- [ ] Specify where to find these classifiers


## Workflow
Run `workflow.rb` to get all candidate sentences to contain a TRI (Transcription Regulation Interation) along their MoR (Mode of Regulation). To run, use:

```
rbbt workflow.rb
# TODO - Miguel - What was the code exactly? 
```

## Postprocessing
Check `scripts/postprocessing/prepare_ExTRI2_resource.ipynb` for an explanation on how was the final ExTRI2 resource created.

## Analysis
`analysis/ExTRI2_results_analysis.ipynb` contains all analysis and figures created for the paper.


# Environments
How to set up the .general_env (used to run all scripts but classifiers training)
```
python3 -m venv .general_env
source .general_venv/bin/activate/
pip install ipykernel
python3 -m ipykernel install --user --name .general_env
pip install pandas matplotlib torch biopython
```

- [ ] Ensure ^ is complete & explain classifiers_training env too
