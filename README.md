# ExTRI2

## Overview
This repo contains the scripts and datasets used for the development of the ExTRI2 pipeline.

## Setup
<span style='color: red'>TODO - Explain how to setup RBBT</span>

## Classifiers training
* Classifiers were trained using the folder `classifiers_training/`. See the `README` inside for an explanation on how to run them.
* The best-performing model was chosen with the notebook `analysis/classifiers_comparison.ipynb`
* The classifiers were used to improve the training dataset through selecting sentences to re-validate that were potentially incorrect. Inside the folder `scripts/classifiers_training/`:
    * See `prepare_reannotation_Excels.ipynb` to see how the reannotation excels were created.
    * See `update_tri_sentences.ipynb` and `make_train_data.ipynb` for the notebooks used to update the dataset & prepare the tsvs for training the models.
* See `classifiers_training/README.md` for an explanation on how models were trained.

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
<span style='color: red'>TODO - specify where to find these classifiers</span>


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

<span style='color: red'>TODO - Ensure ^ is complete & explain classifiers_training env too</span>
