import matplotlib.pyplot as plt
import numpy as np
import glob

def plot_hist(path):   
    dataFile = open(path, 'r')
    # expecting per line: viterbipath;probability
    lines = dataFile.readlines()
    probabilities = [line.split(';')[1] for line in lines]
    probabilities = [float(probabilities[i][:len(probabilities[i])]) for i in range(len(probabilities))]

    n, bins, patches = plt.hist(probabilities, density=True)

    plt.title("Histogram of Viterbi-probabilities")
    plt.xlabel('Viterbi-Probability')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.show()
    
    
def plot_roc(path):
    dataFile = open(path, 'r')
    # expecting: first line: true positive rate
                #second line: false positive rate
    lines = dataFile.readlines()

    TPR = lines[0][1:len(lines[0])-2].split(',')
    TPR = [float(x) for x in TPR]
    FPR = lines[1][1:len(lines[1])-2].split(',')
    FPR = [float(x) for x in FPR]

    auc = np.trapz(TPR, FPR)

    plt.title("ROC-Curve with AUC value:" + str(auc))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.plot(FPR, TPR)
    plt.axis([0.0, 1.0, 0.0, 1.0])
    plt.show() 