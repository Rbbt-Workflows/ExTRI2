import os
import json
from tqdm.auto import tqdm
import torch
import torch.nn as nn
from torch.optim import AdamW

from transformers import AutoModelForSequenceClassification, get_linear_schedule_with_warmup
import lightning as L
from torchmetrics.functional import accuracy, f1_score, auroc, recall, precision, confusion_matrix
from .config import params


class Classifier(L.LightningModule):
    def __init__(self, tokenizer, p: params, train_length: int = 0):
        super().__init__()

        # Define model
        self.model = AutoModelForSequenceClassification.from_pretrained(
            p.pretrained_model_path, return_dict=True,
            num_labels=p.num_labels, id2label=p.id2label, label2id=p.label2id
        )
        # Resize the model's token embeddings (if required)
        if len(p.added_tokens) != 0: 
            self.model.resize_token_embeddings(len(tokenizer))

        self.batch_size = p.batch_size
        self.num_epochs = p.num_epochs
        self.learning_rate = p.learning_rate
        self.train_length = train_length  # Length of training set for the optimizer

        # Define criterion
        if p.metrics_task == 'binary':
            self.criterion = nn.CrossEntropyLoss(weight=torch.tensor(p.train_loss_weights))
        else:
            self.criterion = nn.CrossEntropyLoss()

    def forward(self, input_ids, attention_mask, labels=None):
        output = self.model(input_ids, attention_mask=attention_mask)
        loss = self.criterion(output.logits, labels) if labels is not None else 0
        return loss, output

    def _step(self, batch, batch_idx, loss_log=''):
        input_ids = batch["input_ids"]
        attention_mask = batch["attention_mask"]
        labels = batch["labels"]
        loss, outputs = self(input_ids, attention_mask, labels)
        self.log(loss_log, loss, prog_bar=True, logger=True)
        return loss, outputs, labels
        
    def training_step(self, batch, batch_idx):
        loss, outputs, labels = self._step(batch, batch_idx, loss_log = 'train_loss')
        return {"loss": loss, "predictions": outputs, "labels": labels}

    def validation_step(self, batch, batch_idx):
        loss, outputs, labels = self._step(batch, batch_idx, loss_log = 'val_loss')
        return loss

    def test_step(self, batch, batch_idx):
        loss, outputs, labels = self._step(batch, batch_idx, loss_log = 'test_loss')
        return loss

    def configure_optimizers(self):

        # Calculate total training steps and warmup steps
        steps_per_epoch = self.train_length // self.batch_size
        total_training_steps = steps_per_epoch * self.num_epochs
        n_warmup_steps = total_training_steps // 5

        # Define optimizer and scheduler
        optimizer = AdamW(self.parameters(), lr=self.learning_rate)

        scheduler = get_linear_schedule_with_warmup(
            optimizer,
            num_warmup_steps=n_warmup_steps,
            num_training_steps=total_training_steps,
        )

        return dict(
            optimizer=optimizer, lr_scheduler=dict(scheduler=scheduler, interval="step")
        )

    def calculate_training_steps(train_len, batch_size, num_epochs):
        """
        Calculate the total training steps and warmup steps for a training configuration.
        """
        steps_per_epoch = train_len // batch_size
        total_training_steps = steps_per_epoch * num_epochs
        warmup_steps = total_training_steps // 5
        return warmup_steps, total_training_steps

def evaluate_one_epoch(model, p, k: int, val_dataloader, downsampled=False):
    """
    Evaluate one epoch and add a new row to the report_df
    """

    model.eval()
    model.freeze()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    # Get logits from the model
    logits, labels, texts = [], [], []
    for batch in tqdm(val_dataloader):
        _, output = model(batch["input_ids"].to(device), 
                          batch["attention_mask"].to(device))
        logits.append(output.logits)
        labels.append(batch["labels"])
        texts.extend(batch["texts"])
    
    logits = torch.cat(logits).detach().cpu()
    labels = torch.cat(labels).detach().cpu()

    # Apply softmax to squash logits into [0, 1] range
    probabilities = torch.nn.Softmax(dim=-1)(logits)
    if p.metrics_task == 'binary':
        probabilities = probabilities[:, 1]

    # Get predictions with argmax
    predictions = torch.argmax(logits, dim=-1)
   
    # Compute metrics
    kwargs = {"task": p.metrics_task, "average": p.metrics_average, "num_classes": p.num_labels}
    metrics = {
        'accuracy':  accuracy(  predictions, labels, **kwargs),
        'f1_score':  f1_score(  predictions, labels, **kwargs),
        'recall':    recall(    predictions, labels, **kwargs),
        'precision': precision( predictions, labels, **kwargs),
        'auroc':     auroc(   probabilities, labels, **kwargs),
    }

    # Add entries to p.report_df of params
    if p.metrics_task == 'binary':
        for metric_name, metric_val in metrics.items():
            p.report_df[metric_name].append(round(metric_val.item(), 5))
    else:
        for metric_name, metric_val in metrics.items():
            p.report_df[metric_name].append(metric_val.tolist())

    p.report_df['fold'].append(k+1)
    p.report_df['confusion'].append(confusion_matrix(probabilities, labels, task=p.metrics_task, num_classes = p.num_labels).tolist())
    
    # Save probabilities and labels for further analysis
    p.report_df['probabilities'].append(probabilities.tolist())
    p.report_df['labels'].append(labels.tolist())
    p.report_df['val_texts'].append(texts)
    
    return

def save_model(p, model, tokenizer):
    '''
    Save model, tokenizer and params (as train_params.json) in p.model_path 
    '''
    if not os.path.exists(p.model_path):
        os.mkdir(p.model_path)
    
    model.model.save_pretrained(p.model_path)
    tokenizer.save_pretrained(p.model_path)
    
    with open(os.path.join(p.model_path, "train_params.json"), "w") as outfile:
        outfile.write(json.dumps(vars(p), indent=4))
    
    print("Model has been saved")

    return