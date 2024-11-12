# Train the MoR (Mode of Regulation) classifier
# Run on the HPC by calling it through 'run_classifiers.sh'
# To run independently: python MoR_train.py [--kwargs]
# Use -h to get a list of the parameters


from src.train import MoR_train

if __name__ == '__main__':
    MoR_train()



    
