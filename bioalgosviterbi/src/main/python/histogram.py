import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        print("Please provide data file")
    
    dataFile = open(sys.argv[1], 'r')
    # expecting per line: viterbipath;probability
    lines = dataFile.readlines()
    probabilities = [line.split(';')[1] for line in lines]
    probabilities = [float(probabilities[i][:len(probabilities[i])]) for i in range(len(probabilities))]

    n, bins, patches = plt.hist(probabilities, density=True, bins=5)

    plt.title("Histogram of Viterbi-probabilities")
    plt.ylabel('Viterbi-Probability')
    plt.xlabel('Probability')
    plt.grid(True)
    plt.show()