import argparse

class params:
    def __init__(self):
        self.class_name         = self.__class__.__name__
        self.description        = ''
        self.models_path        = ''
        self.data_folder        = ''
        self.pretrained_model_folder = 'pretrained_models/'
        self.pretrained_model   = 'distilbert'  # Check 'pretrained_model_path' property to check the options
        self.debug              = False
        self.num_epochs         = 3
        self.learning_rate      = 3e-5
        self.batch_size         = 8
        self.trunc_max_length   = 512
        self.dl_num_workers     = 15
        self.train_loss_weights = [0.5, 0.5]
        self.fold_validation    = True          # Whether to do a 5-fold cross-validation on the train-val split (85% of the data). Else, only train on 1 fold.
        self.get_final_model    = False         # Train on 85% of data and validate with the remaining 15%.
        self.improve_dataset    = False         # Do K-fold validation with 100% of the dataset (no test dataset) to get predictions for all the dataset
        self.metrics_average    = 'micro'  
        self.metrics_task       = 'binary'  
        self.num_labels         = 2
        self.labels_list        = None
        self.added_tokens       = []
        self.report_df          = {key: [] for key in ['fold', 'accuracy', 'f1_score', 
                                                       'recall', 'precision', 'auroc', 
                                                       'confusion', 'probabilities', 'labels', 'val_texts']}

    @property
    def model_path(self):
        return self.models_path + self.model_name
    
    @property
    def pretrained_model_path(self):
        full_name = self.pretrained_model_folder    
        # If not, use a pretrained model from transformers
        if self.pretrained_model == 'distilbert':
            full_name += 'distilbert-base-uncased'
        elif self.pretrained_model == 'BiomedNLP':
            full_name += 'BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext'
        elif self.pretrained_model == 'BiomedNLP-large':
            full_name += 'BiomedNLP-BiomedBERT-large-uncased-abstract'
        elif self.pretrained_model == 'BioLinkBERT':
            full_name += 'BioLinkBERT-base'
        elif self.pretrained_model == 'BioBERT':
            full_name += 'biobert-base-cased-v1.2'
        elif self.pretrained_model == 'sciBERT':
            full_name += 'scibert_scivocab_uncased'
        else: 
            raise ValueError(f"Invalid pretrained model name '{self.pretrained_model}'. Valid names are 'distilbert' and 'BiomedNLP'.")

        return full_name
   
    @property
    def pretrained_model_abbreviation(self):
        if self.model_name == '':
            abbreviation = 'distilbert'
        elif self.model_name == '':
            abbreviation = ''
        else:
            raise ValueError
        return abbreviation

def create_parser_for_params(params):
    parser = argparse.ArgumentParser(description='Train models with dynamic parameters.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Decide the type of training to do

    # List of parameters to exclude from command-line arguments
    exclude_params = ['class_name', 'train_type', 'models_path', 'dataset_path', 'id2label', 'label2id', 'num_labels', 'labels_list', 'epochs_metrics']

    for attribute, value in vars(params).items():
        if attribute.startswith('_') or attribute in exclude_params:
            continue  # Skip private attributes and excluded parameters

        if value is None:
            parser.add_argument(f'--{attribute}', default=value, help=f'Type inferred. Value: {value}')
        elif isinstance(value, list):
            parser.add_argument(f'--{attribute}', nargs='+', default=value)
        elif isinstance(value, bool):
            parser.add_argument(f'--{attribute}', type=lambda x: (str(x).lower() in ['true', '1', 'yes']), default=value, help=f'Boolean. Value: {value}')
        else:
            parser.add_argument(f'--{attribute}', type=type(value), default=value, help=f'Type: {type(value).__name__}. Value: {value}')

    args = parser.parse_args()
    return args

def update_params(class_type, args) -> params:
    '''Create an instance of class_type and return it with updated parameters'''
    
    p = class_type()
    
    # Update the params instance with arguments provided, ensuring exact match
    for arg in vars(args):
        if hasattr(p, arg):
            setattr(p, arg, getattr(args, arg))

    return p

class TRI_params(params):
    def __init__(self):
        super().__init__()
        self.models_path        = 'saved_models/TRI_classifier/'
        self.data_folder        = 'data/'
        self.trunc_max_length   = 512
        self.metrics_average    = 'binary'        
        self.id2label           = {0: "DISCARD", 1: "ACCEPT"}
        self.label2id           = {"ACCEPT": 0, "DISCARD": 1}
        self.num_labels         = 2
        self.masked             = False
        self.span               = False

    @property
    def model_name(self):
        
        model = self.pretrained_model
        span = '_span' if self.span else ''
        final = '_final' if self.get_final_model else ''
        tokens = '_' + '_'.join(self.added_tokens) if len(self.added_tokens) != 0 else ''
        fold = '_1fold' if not self.fold_validation else ''
        epochs = f'_{self.num_epochs}e' if self.num_epochs != 3 else ''
        debug = '_debug' if self.debug else ''
        learning_rate = f'_lr{self.learning_rate}'

        return f'{model}{span}{tokens}{fold}{final}{epochs}{learning_rate}{debug}'
    
class MoR_params(params):
    def __init__(self):
        super().__init__()
        self.models_path        = 'saved_models/MoR_classifier/'
        self.data_folder       = 'data/'
        self.added_tokens       = []
        self.trunc_max_length   = 512
        self.id2label           = {0: "UNDEFINED", 1: "ACTIVATION", 2: "REPRESSION"}
        self.label2id           = {"UNDEFINED": 0, "ACTIVATION": 1, "REPRESSION": 2}
        self.num_labels         = 3
        self.metrics_average    = 'none'  
        self.metrics_task       = 'multiclass'
        self.labels_list        = [0,1,2]
        self.span               = False

    @property
    def model_name(self):
        
        model = self.pretrained_model
        tokens = '_' + '_'.join(self.added_tokens) if len(self.added_tokens) != 0 else ''
        span = '_span' if self.span else ''
        fold = '_1fold' if not self.fold_validation else ''
        final = '_final' if self.get_final_model else ''
        epochs = f'_{self.num_epochs}e' if self.num_epochs != 3 else ''        
        debug = '_debug' if self.debug else ''
        learning_rate = f'_lr{self.learning_rate}'

        return f'{model}{span}{tokens}{fold}{final}{epochs}{learning_rate}{debug}'
