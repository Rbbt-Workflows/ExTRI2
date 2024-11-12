import pandas as pd
from torch.utils.data import Dataset, DataLoader
import lightning as L
from sklearn.model_selection import KFold, train_test_split
from .config import params
import sys

class GenericClassifierDataset(Dataset):
    def __init__(
        self, data: pd.DataFrame, tokenizer, trunc_max_length=512
    ):
        self.tokenizer = tokenizer
        self.data = data
        self.trunc_max_len = trunc_max_length

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index: int):
        data_row = self.data.iloc[index]
        texts = data_row.texts
        labels = data_row.labels
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
            texts=texts,
            labels=labels,
            input_ids=encoding["input_ids"].flatten(),
            attention_mask=encoding["attention_mask"].flatten(),
        )

class KFoldDataModule(L.LightningDataModule):
    def __init__(self, train_val_df, test_df, tokenizer, k: None, p: params):
        super().__init__()
        self.p = p
        self.k = k      # Fold to train with 
        self.batch_size = p.batch_size
        self.train_val_df = train_val_df
        self.test_df = test_df
        self.tokenizer = tokenizer
        self.trunc_max_len = p.trunc_max_length
        self.dl_num_workers = p.dl_num_workers
        self.params_class = p.class_name  # TRI_params or MoR_params
        self.added_tokens = p.added_tokens

    def setup(self, stage=None):

        if self.p.get_final_model:
            # Validate with the test set
            train_df, val_df = self.train_val_df, self.test_df

        else:
            # Create 5 folds from the train_val_df. Train with the fold number k.
            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            all_splits = [k for k in kf.split(self.train_val_df)]
            train_idx, val_idx = all_splits[self.k]
            train_idx, val_idx = train_idx.tolist(), val_idx.tolist()
            train_df, val_df = (
                self.train_val_df.iloc[train_idx],
                self.train_val_df.iloc[val_idx],
            )
                
        self.train_dataset = GenericClassifierDataset(train_df, self.tokenizer, self.trunc_max_len)
        self.val_dataset =   GenericClassifierDataset(val_df, self.tokenizer, self.trunc_max_len)
        self.test_dataset =  GenericClassifierDataset(self.test_df, self.tokenizer, self.trunc_max_len)

    def _create_dataloader(self, dataset, shuffle):
        """Common dataloader for train, test and val"""
        return DataLoader(
            dataset,
            batch_size=self.batch_size,
            shuffle=shuffle,
            num_workers=self.dl_num_workers,
        )
    
    def train_dataloader(self):
        return self._create_dataloader(self.train_dataset, shuffle=True)

    def val_dataloader(self):
        return self._create_dataloader(self.val_dataset, shuffle=False)

    def test_dataloader(self):
        return self._create_dataloader(self.test_dataset, shuffle=False)

def load_datasets(p):
    if p.class_name == 'TRI_params':
        if p.masked is True:
            data = pd.read_csv(p.data_folder + 'TRI_masked.tsv',sep='\t')
        elif p.span is True:
            data = pd.read_csv(p.data_folder + 'TRI_span_data.tsv',sep='\t')
        else:
            data = pd.read_csv(p.data_folder + 'TRI_data.tsv',sep='\t')
        train_data, test_data = train_test_split(data, test_size=0.15, random_state=42)

    elif p.class_name == 'MoR_params':
        if p.span is True:
            data = pd.read_csv(p.data_folder + 'MoR_span_data.tsv',sep='\t')
        else:
            data = pd.read_csv(p.data_folder + 'MoR_data.tsv',sep='\t')
        train_data, test_data = train_test_split(data, test_size=0.15, random_state=42)
        
    else:
        raise ValueError(f'Incorrect class name: {p.class_name}')

    if p.improve_dataset:
        # Use the whole dataset for training, to get all incorrect rows.
        train_data = data

    if p.debug == True:
        train_data, _ = train_test_split(train_data, test_size=0.99, random_state=42)
        test_data, _  = train_test_split(test_data, test_size=0.99, random_state=42)


    return train_data, test_data