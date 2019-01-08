package viterbi;

public class Viterbi {
    private static double[][] toLog(double[][] matrix) {
        for(int i = 0; i<matrix.length; i++) {
            for(int j = 0; j<matrix[0].length; j++) {
                matrix[i][j] = Math.log(matrix[i][j]);
            }
        }
        return matrix;
    }



    public static int[] calc(int[] observations, double[] initialProbabilities, double[][] transitionMatrix, double[][] emissionMatrix) {
        int countStateSpace = transitionMatrix.length;
        double[][] viterbiVar = new double[countStateSpace][observations.length];
        int[][] backtrackingVar = new int[countStateSpace][observations.length];

        double[][] logTransitionMatrix = Viterbi.toLog(transitionMatrix.clone());
        double[][] logEmissionMatrix = Viterbi.toLog(emissionMatrix.clone());

        for(int state = 0; state < countStateSpace; state++) {
            viterbiVar[state][0] = Math.log(initialProbabilities[state]) + logEmissionMatrix[state][observations[0]];
        }

        for(int observation = 1; observation < observations.length; observation++) {
            for(int state = 0; state < countStateSpace; state++) {
                double maxScore = viterbiVar[0][observation-1] + logTransitionMatrix[0][state];
                int argMaxScore = 0;
                for(int i = 1; i < countStateSpace; i++) {
                    double score = viterbiVar[i][observation-1] + logTransitionMatrix[i][state];
                    if (score > maxScore) {
                        maxScore = score;
                        argMaxScore = i;
                    }
                }
                maxScore += logEmissionMatrix[state][observations[observation]];
                viterbiVar[state][observation] = maxScore;
                backtrackingVar[state][observation] = argMaxScore;
            }
        }
        int[] path = new int[observations.length];
        int last = path.length - 1;
        path[last] = backtrackingVar[0][last];
        for(int state = 1; state<countStateSpace; state++) {
            if (backtrackingVar[state][last] > path[last]) {
                path[last] = backtrackingVar[state][last];
            }
        }
        for(int observation = last; observation > 0; observation--) {
            path[observation-1] = backtrackingVar[path[observation]][observation];
        }
        return path;
    }
}
