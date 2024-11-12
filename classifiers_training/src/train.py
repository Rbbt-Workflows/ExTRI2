import torch
import pandas as pd
import lightning as L
from transformers import AutoTokenizer
from lightning.pytorch.callbacks import ModelCheckpoint, EarlyStopping
from lightning.pytorch.loggers import TensorBoardLogger
from sklearn.model_selection import KFold
from .config import params, TRI_params, MoR_params, create_parser_for_params, update_params
from .dataset import KFoldDataModule, load_datasets
from .model import Classifier, evaluate_one_epoch, save_model

RANDOM_SEED = 42
L.seed_everything(RANDOM_SEED)
torch.set_float32_matmul_precision('medium') # to make pytorch happy. 'high' will give a higher precision but slower running

def train(p):

    # Load datasets
    train_data, test_data = load_datasets(p)

    # Define the kfold splits
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    all_splits = [k for k in kf.split(train_data)]

    # Use only 1 fold if getting the final model or p.fold_validation = False
    all_splits = all_splits[:1] if (p.get_final_model or (not p.fold_validation)) else all_splits
 
    # Define the tokenizer
    tokenizer = AutoTokenizer.from_pretrained(p.pretrained_model_path)
    
    # Modify tokenizer (if required)
    if len(p.added_tokens) != 0:
        new_num_tokens = tokenizer.add_tokens(p.added_tokens)
        print(f'{new_num_tokens} tokens have been added to the tokenizer: {p.added_tokens}')

    # Train and validate for each split
    for k, split in enumerate(all_splits):

        # Checkpoints
        checkpoint_callback = ModelCheckpoint(
            dirpath=p.model_path + '/checkpoints',
            filename=f"best-checkpoint_{k}",
            save_top_k=1,
            verbose=True,
            monitor="val_loss",
            mode="min",
        )

        # Early stopping triggers
        early_stopping_callback = EarlyStopping(monitor="val_loss", patience=2)
        
        # Define the trainer
        trainer = L.Trainer(
            callbacks=[early_stopping_callback, checkpoint_callback],
            max_epochs=p.num_epochs,
            default_root_dir=p.model_path,
            # fast_dev_run=True,  # Run only 1 batch per epoch to debug (doesn't save the results)
        )

        # Define data module for the fold k
        data_module = KFoldDataModule(train_data, test_data, tokenizer, k=k, p=p)

        # Build the model
        train_length = len(train_data) if p.get_final_model else len(split[0])
        model = Classifier(tokenizer, p=p, train_length=train_length)
        # Train
        trainer.fit(model=model, datamodule=data_module)

        # Evaluate
        evaluate_one_epoch(model, p, k, data_module.val_dataloader())

    report_df = pd.DataFrame(p.report_df)
    report_df.set_index("fold", inplace=True)

    # Print results
    if 'val_type' in report_df.columns:
        print(report_df[['val_type', 'accuracy', 'f1_score', 'recall', 'precision', 'auroc']])
    else:
        print(report_df[['accuracy', 'f1_score', 'recall', 'precision', 'auroc']])


    # Save model, tokenizer and params
    save_model(p, model, tokenizer)

    return    

def TRI_train():
    # Parse the command-line arguments
    args = create_parser_for_params(TRI_params())

    # Update the params instance with the arguments
    p = update_params(TRI_params, args)

    # Train and save
    train(p) 

def MoR_train():
    # Parse the command-line arguments
    args = create_parser_for_params(MoR_params())

    # Update the params instance with the arguments
    p = update_params(MoR_params, args)

    # Train and save
    train(p) 

if __name__ == "__main__":
    train(params())
