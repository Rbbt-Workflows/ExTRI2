'''
Generates predictions for a given dataset using the specified model. Appends the predictions as a new column to the dataset.
Requires: paths to the dataset and the model, class_type (TRI, MoR) and output_file.
Dataset must contain a column named 'texts' with the texts to classify

Designed for standalone use, without dependencies on internal code.

Example use: 
python src/predict.py --class_type TRI --model_path final_models/TRI_model/ --data_path <input.tsv> --output_path <output.tsv>
'''

import torch
from torch.utils.data import Dataset, DataLoader
from tqdm.auto import tqdm
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import pandas as pd
import argparse

def setup_parser():
    """Configure and return the argument parser."""
    parser = argparse.ArgumentParser(description='Predict using a trained model.')
    parser.add_argument('--class_type', type=str, required=True, 
                        choices=['TRI', 'MoR'], 
                        help='Class type must be one of: TRI, MoR.')
    parser.add_argument('--model_path', type=str, required=True, 
                        help='Path to the trained model file.')
    parser.add_argument('--data_path', type=str, required=True, 
                        help='Path to the input data for prediction.')
    parser.add_argument('--output_path', type=str, required=True, 
                        help='Path to save the predictions.')

    return parser

def load_model_tokenizer_from_saved(model_path):
    model = AutoModelForSequenceClassification.from_pretrained(model_path)
    tokenizer = AutoTokenizer.from_pretrained(model_path)    
    return model, tokenizer

class PredictDataset(Dataset):
    
    def __init__(
        self, data: pd.DataFrame, tokenizer: AutoTokenizer, trunc_max_length=512
    ):
        self.tokenizer = tokenizer
        self.data = data
        self.trunc_max_len = trunc_max_length

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index: int):
        data_row = self.data.iloc[index]
        texts = data_row.texts
        encoding = self.tokenizer.encode_plus(
            texts,
            add_special_tokens=True,
            max_length=self.trunc_max_len,
            return_token_type_ids=False,
            padding="max_length",
            truncation=True,
            return_attention_mask=True,
            return_tensors="pt",
        )
        return dict(
            input_ids=encoding["input_ids"].flatten(),
            attention_mask=encoding["attention_mask"].flatten(),
        )

def predict(dataset: pd.DataFrame, model_path: str, id2label: dict):
    '''
    Generates predictions for a given dataset using a pretrained model, and appends these predictions as a new column to the input DataFrame.

    Parameters:
    - dataset (pd.DataFrame): Must contain the column 'texts'
    - model_path (str): Used to load the model and tokenizer.
    - id2label (dict): A dictionary mapping the model's output indices to their corresponding labels.

    Returns:
    - None: The function modifies the input DataFrame in-place by adding a 'predictions' column with the labeled predictions.
    '''

    # Get model, tokenizer and dataloader
    model, tokenizer = load_model_tokenizer_from_saved(model_path)
    dataloader = DataLoader(
        PredictDataset(dataset, tokenizer = tokenizer),
        batch_size=16,
        shuffle=False,
    )

    # Put the model in evaluation mode
    model.eval()

    # Get the predictions
    predictions_list = []
    probabilities_list = []
    
    progress_bar = tqdm(total=len(dataloader))
    for batch in dataloader:
        with torch.no_grad():
            output = model(**batch)
        logits = output.logits
        preds = torch.argmax(logits, dim=-1)
        probs = torch.nn.Softmax(dim=-1)(logits)
        predictions_list.extend(preds.tolist())  # Assuming outputs is a tensor
        probabilities_list.extend([p[1] for p in probs.tolist()])
        progress_bar.update(1)

    dataset['predictions'] = predictions_list
    dataset['predictions'] = dataset['predictions'].replace(id2label)
    dataset['probabilities'] = probabilities_list

    return

def main():
    # Get CLI arguments
    parser = setup_parser()
    args = parser.parse_args()

    # Mapping from class_type to id2label
    id2label_map = {
        'TRI':  {0: "Not valid", 1: "Valid"},
        'MoR': {0: "UNDEFINED", 1: "ACTIVATION", 2: "REPRESSION"}
    }
    id2label = id2label_map[args.class_type]

    # Create dataset
    dataset = pd.read_csv(args.data_path, sep='\t', header=0)

    # Predict and save results as a new column in the dataset
    predict(dataset, args.model_path, id2label)

    # Save the dataset in the output file
    dataset.to_csv(args.output_path, sep='\t')

if __name__ == "__main__":
    main()