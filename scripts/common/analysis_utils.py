# Functions used for plotting and visualising the data

import matplotlib.pyplot as plt
import pandas as pd

def prettify_plots(fontsize=12, titlesize=12, labelsize=10, linewidth=1):
    '''
    prettify the plots. Set multiple rcParams
    '''
    plt.rcParams.update({
    'font.size': fontsize,                # Default font size
    'axes.titlesize': titlesize,          # Font size of the plot titles
    'axes.labelsize': labelsize,          # Font size of the axis labels
    'xtick.labelsize': labelsize,         # Font size of the x-axis tick labels
    'ytick.labelsize': labelsize,         # Font size of the y-axis tick labels
    'legend.fontsize': labelsize,         # Font size of the legend
    'lines.linewidth': linewidth,         # Default line width
    })
    return


# 終 TODO　- this plotting is horrible. I should probably remove it and create more useful ones.
#            try moving it to seaborn?
def plot_seqlen_distribution(texts, labels, suptitle, red_title="Discarded Data", green_title="Accepted Data"):
    """
    Plot the distribution of sequence lengths in a dataset, separated by labels 0 and 1.

    Parameters:
    - texts (list): A list of text data.
    - labels (list): A list of labels corresponding to the text data.
    - suptitle (str): The main title of the plot.
    - red_title (str, optional): The title for data with label=0 (default: "Discarded Data").
    - green_title (str, optional): The title for data with label=1 (default: "Accepted Data").

    Returns:
    - None: The function generates and displays the plot.
    """
    data = pd.DataFrame()
    data["Text"] = texts
    data["Label"] = labels
    
    # Plot the data distribution
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    plt.suptitle(suptitle)
    
    # Plot the distribution of input ids
    for i, color, title in [[0, 'red', red_title], [1, 'green', green_title]]:
        list = data[data['Label'] == i ]['Text'].apply(len)
        axs[i].hist(list, bins=20, edgecolor='k', color=color, alpha=0.7)
        axs[i].set_xlabel('Sequence Length')
        axs[i].set_ylabel('Frequency')
        suptitle = f"{len(data[data['Label'] == i])} abstracts ({len(data[data['Label'] == i])/len(data)*100:.2f}%)"
        axs[i].set_title(f"{title}\n{suptitle}")  
    axs[0].hist(data[data['Label'] == 1 ]['Text'].apply(len), bins=20, edgecolor='k', color="green", alpha=0.1)
    
    plt.tight_layout()
    return

def find_false_predictions(val_results_dict: dict, predictions: list, probabilities: list, true_labels: list, texts: list):
    '''
    Find and format false positive and false negative examples with their probabilities.
    '''

    # Get indices wrongly predicted texts
    # False positive: pred=1 but label=0. pred>label. Opposite for false negatives
    false_pos_idx = [i for i,(pred,label) in enumerate(zip(predictions, true_labels)) if pred > label ]
    false_neg_idx = [i for i,(pred,label) in enumerate(zip(predictions, true_labels)) if pred < label ]
    
    # Get the wrongly predicted texts and their probabilities. Sort them
    false_pos_texts = [[texts[i], f"{probabilities[i][1]:.4f}"] for i in false_pos_idx]
    false_neg_texts = [[texts[i], f"{probabilities[i][1]:.4f}"] for i in false_neg_idx]

    # Sort them
    false_pos_texts = sorted(false_pos_texts, key=lambda x: x[1], reverse=True)
    false_neg_texts = sorted(false_neg_texts, key=lambda x: x[1])
    
    # Add false positives and false negatives to the dictionary
    val_results_dict['false_positives'] = false_pos_texts
    val_results_dict['false_negatives'] = false_neg_texts

    return false_pos_texts, false_neg_texts