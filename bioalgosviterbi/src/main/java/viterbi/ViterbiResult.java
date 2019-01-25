package viterbi;

import java.util.Arrays;

public class ViterbiResult {
    private final int[] viterbiPath;
    private final double maxProbability;

    public ViterbiResult(int[] viterbiPath, double maxProbability) {
        this.viterbiPath = viterbiPath;
        this.maxProbability = maxProbability;
    }

    @Override
    public String toString(){
        return Arrays.toString(this.viterbiPath) + ";" + maxProbability;
    }

    public int[] getViterbiPath() {
        return this.viterbiPath;
    }

    public double getMaxProbability() {
        return this.maxProbability;
    }
}