# ExTRI2

## Overview
ExTRI2 identifies candidate sentences that describe TF–target interactions from PubTator3-annotated abstracts, classifies whether they are valid TRIs, and assigns an MoR label (Activation, Repression, Undefined).

This repo contains the ExTRI2 workflow, training code for the TRI/MoR classifiers, preprocessing/postprocessing utilities, and analysis used in the paper.

## Setup
Create a Python environment and install dependencies for the main pipeline. We used python=3.12

```
python -m venv .extri2_venv
source .extri2_venv/bin/activate
pip install -r requirements.txt
python -m ipykernel install --user --name=extri2_venv
```

RBBT setup: please follow the official installation docs for your platform. Once installed, you can run workflows from this repository. Useful starting points:
- RBBT: https://github.com/mikisvaz/rbbt
- Workflows guide: https://github.com/mikisvaz/rbbt/wiki/Workflows


## Quickstart
1) **Prepare inputs for ExTRI2**
```
# Build TF list tf_entrez_code.list
cd scripts/preprocessing/
jupyter nbconvert --to notebook --execute get_NCBI_TF_IDs.ipynb

# Download PubTator3 files and prepare them for ExTRI2
cd ../../
./scripts/preprocessing/get_all_pubtators.sh
python3 scripts/preprocessing/prepare_pubtator_for_ExTRI2.py
```

2) **Provide models**
Train models following `classifiers_training/README.md`, or use existing ones. ExTRI2 can accept either:
- Model names registered in RBBT (e.g. `TRI_model`, `MoR_model`), or
- Absolute paths to the exported model directories.

3) **Run the workflow**
List available tasks and options:
```
rbbt workflow.rb --help
```

Run on a single PubTator file (example):
```
rbbt workflow.rb tri_candidates --pubtator_file examples/tri_candidates/pubtator_file
rbbt workflow.rb tri_sentences --pubtator_file examples/tri_candidates/pubtator_file --tri_model /path/to/TRI_model
rbbt workflow.rb tri_MoR --pubtator_file examples/tri_candidates/pubtator_file --mor_model /path/to/MoR_model
```

Batch run over all `.pubtator` files detected by your RBBT data path:
```
rbbt workflow.rb ExTRI2
```

4) **Postprocessing**
Convert workflow outputs into the final ExTRI2 dataset and apply renormalisations:
```
cd scripts/postprocessing/
python postprocessing.py

cd ../../analysis/
jupyter nbconvert --to notebook --execute repo_to_paper.ipynb

```

## Training the Classifiers
- Environment (Python 3.11 recommended for training):

```
cd classifiers_training/
python3 -m venv .clf_venv
source .clf_venv/bin/activate
pip install -r requirements.txt
python3 -m ipykernel install --user --name extri2-classifiers
```

- Train and compare models:

```
python TRI_train.py -h
python MoR_train.py -h
```

See `analysis/classifiers_comparison.ipynb` for model selection and `scripts/classifiers_training/` notebooks to prepare and refine training data:
    - `prepare_reannotation_Excels.ipynb` prepares the sentences to revise.
    - `update_tri_sentences.ipynb` and `make_train_data.ipynb` update the datasets and prepare the files used for training the models.
- See `classifiers_training/README.md` for details on training and model selection.


## Steps to obtain the results:
The final ExTRI2 dataset required training the classifiers and improving the training dataset, preparing files for the `workflow.rb`, postprocessing the resulting files, and preparing sentences for validation. This was achieved by running the following scripts:


## Analysis
`analysis/repo_to_paper.ipynb` contains all analysis and figures created for the paper, as well as links to scripts used, divided by each of the paper's sections.


## Repository structure
- `workflow.rb`: RBBT workflow to obtain TRI sentences and MoR labels from PubTator files.
- `classifiers_training/`: Code and data to train the TRI and MoR classifiers. See `classifiers_training/README.md` for details.
- `scripts/`: Helper scripts.
  - `preprocessing/`: Obtain TF IDs, fetch PubTator files, and prepare input.
  - `postprocessing/`: Convert workflow outputs into the final ExTRI2 dataset and perform renormalisations.
  - `validation/`: Prepare and summarise manual validation sets.
- `data/`: Raw and intermediate files required to run the workflow (see `data/README.md`).
- `analysis/`: Notebooks and scripts for figures and comparisons used in the paper.


## Citation
If you use ExTRI2 in your work, please cite:
Fàbrega, N., Thommesen, L., Valencia, A., Lægreid, A., Vazquez, M.  
*ExTRI2: Harnessing Transformers for Enhanced Mining of Transcription Regulation Interactions* (under review).

## License
This project is licensed under the terms in `LICENSE`.