# Train the TRI (Transcription Regulation Interaction) classifier
# Run on the HPC by calling it through 'run_classifiers.sh'
# To run independently: python TRI_train.py [--kwargs]
# Use -h to get a list of the parameters


from src.train import TRI_train

if __name__ == '__main__':
    TRI_train()

    
