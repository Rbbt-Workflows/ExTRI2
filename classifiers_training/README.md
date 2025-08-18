# Classifiers training
This folder contains the scripts and data used to train the TRI and MoR classifiers, as well as the parameters and results obtained from training them.

## Structure
Scripts:
* `src/` contains the scripts to train the classifiers. It also includes `predict.py`, a standalone script for the prediction of datasets using the trained models. 
* `MoR_train.py` and `TRI_train.py` are the scripts to run to train the classifiers. 
* `run_classifiers.sh` used to send jobs to run to the MN5 HPC in batch.
Datasets and models:
* `data/` contains the datasets used to train the classifiers.
* `pretrained_models/` contains the pretrained models downloaded from Huggingface (follow `pretrained_models/README.md` to download them)
* `final_models/` contains the final TRI and MoR models used to create the ExTRI2 pipeline.
* `saved_models/` contains a folder for each of the trained models. Inside each, `train_params.json` contains the parameters used & metrics obtained for training
  * The subfolders `before_3rd_iteration`, `before_2nd_iteration` and `before_1st_iteration` contain the models trained with previous versions of the dataset. See  [Dependencies with other folders](#dependencies-with-other-folders) section for more details.


## Create the environment
The environment has been created with python 3.11. It can be installed through:
```
pip install -r requirements.txt
```

## Training the classifiers
Use `TRI_train.py` to train the TRI classifier or `MoR_train.py` to train the MoR classifier. Example:
```
python TRI_train.py --pretrained_model BioLinkBERT --span True --num_epochs 5,3:00:00
```

To get a list of all possible parameters:
```
python TRI_train.py -h
```

Alternatively, to train in batch different models with multiple parameters using the HPC, use:
```
./run_classifiers.sh
```

## Dependencies with other folders
This folder can be run independently from the rest of the folders. 

However:
* The content of `data/` is generated using `../scripts/classifiers_training/make_train_data.ipynb` from the dataset `data/tri_sentences.tsv`, which in turn is created from updating the original dataset `../data/external/original_tri_sentences.tsv` with the reannotated datasets found in `../data/external/dataset_improvement/`.
* The comparison of the different models inside a iteration to find the best-performing one was made in `../analysis/classifiers_comparison.ipynb`. There, the model used in each iteration to find the sentences to reannotate is also specified.
