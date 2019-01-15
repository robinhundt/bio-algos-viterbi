package viterbi;

public class Viterbi {
    /**
     * Calculates the viterbi path given a sequence of observations and a model.
     * The observations are expected to be elements of [0, |EmissionStates|).
     * @param observations Sequence of observations
     * @param initialProbabilities For k element of States, initialProbabilities[k] is the probability that a
     *                             state sequence starts with k
     * @param transitionMatrix Provides the probabilities with which state transitions occur.
     *                         transitionMatrix[j][k] = 0.25 means that there is a chance of 0.25 to go from state
     *                         j to state k. Size must be |States| x |States| .
     * @param emissionMatrix  Provides the probabilities with which an emission occurs, given a state.
     *                        emissionMatrix[j][o] = 0.25 means that there is a 0.25 chance to emit an observation o
     *                        when in state j
     * @return the viterbi path
     */
    public static int[] calc(int[] observations, double[] initialProbabilities, double[][] transitionMatrix, double[][] emissionMatrix) {
        int countStateSpace = transitionMatrix.length;
        double[][] viterbiVar = new double[countStateSpace][observations.length+1];
        int[][] backtrackingVar = new int[countStateSpace][observations.length+1];

        viterbiVar[0][0] = 1;
        
        /*
            Logarithm is applied element wise to the transition and emission matrices. Computing the logithms of
            the probabilities instead of the probabilities themselves ensures that the values stay in ranges which
            can be represented by doubles.
         */
        double[][] logTransitionMatrix = Viterbi.toLog(transitionMatrix.clone());
        double[][] logEmissionMatrix = Viterbi.toLog(emissionMatrix.clone());

        /*
            Initialise first column
         */
        /*
        for(int state = 0; state < countStateSpace; state++) {
            viterbiVar[state][0] = Math.log(initialProbabilities[state]) + logEmissionMatrix[state][observations[0]];
        }
        */

        calcViterbiBacktrackVars(observations, countStateSpace, viterbiVar, backtrackingVar, logTransitionMatrix, logEmissionMatrix);

        return reconstructOptimalPath(observations, countStateSpace, backtrackingVar, viterbiVar);
    }


    private static int[] reconstructOptimalPath(int[] observations, int countStateSpace, int[][] backtrackingVar, double[][] viterbiVar) {
        int[] path = new int[observations.length];
        int last = path.length - 1;
        path[last] = 0;
        double maxLastColumn = viterbiVar[0][last+1];
        for(int state = 1; state<countStateSpace; state++) {
            if (viterbiVar[state][last+1] > maxLastColumn) {
                maxLastColumn = viterbiVar[state][last+1];
                path[last] = state;
            }
        }
        for(int observation = last; observation > 0; observation--) {
            path[observation-1] = backtrackingVar[path[observation]][observation+1];
        }
        return path;
    }

    private static void calcViterbiBacktrackVars(int[] observations, int countStateSpace, double[][] viterbiVar, int[][] backtrackingVar, double[][] logTransitionMatrix, double[][] logEmissionMatrix) {
        for(int observationIdx = 1; observationIdx <= observations.length; observationIdx++) {
            for(int state = 0; state < countStateSpace; state++) {
                double maxScore = viterbiVar[0][observationIdx -1] + logTransitionMatrix[0][state];
                int argMaxScore = 0;
                for(int i = 1; i < countStateSpace; i++) {
                    double score = viterbiVar[i][observationIdx -1] + logTransitionMatrix[i][state];
                    if (score > maxScore) {
                        maxScore = score;
                        argMaxScore = i;
                    }
                }
                maxScore += logEmissionMatrix[state][observations[observationIdx-1]];
                viterbiVar[state][observationIdx] = maxScore;
                backtrackingVar[state][observationIdx] = argMaxScore;
            }
        }
    }

    /**
     * Transforms the provided matrix by an element wise natural logarithm in place.
     * @param matrix 2d matrix that should be transformed.
     * @return the transformed matrix.
     */
    private static double[][] toLog(double[][] matrix) {
        for(int i = 0; i<matrix.length; i++) {
            for(int j = 0; j<matrix[0].length; j++) {
                matrix[i][j] = Math.log(matrix[i][j]);
            }
        }
        return matrix;
    }
}
