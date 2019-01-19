import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        print("Please provide data file")
    
    dataFile = open(sys.argv[1], 'r')
    # expecting: first line: true positive rate
                #second line: false positive rate
    lines = dataFile.readlines()

    TPR = lines[0][1:len(lines[0])-2].split(',')
    TPR = [int(x) for x in TPR]
    FPR = lines[1][1:len(lines[1])-2].split(',')
    FPR = [int(x) for x in FPR]

    auc = np.trapz(TPR, FPR)

    plt.title("ROC-Curve with AUC value:" + str(auc))
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.plot(FPR, TPR)
    plt.axis([0.0, 1.0, 0.0, 1.0])
    plt.show() 