package viterbi;

import phmm.ProfileHMM;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class Viterbi {
    /**
     * Calculates the viterbi path given a sequence of observations and a model.
     * The observations are expected to be elements of [0, |EmissionStates|).
     *
     * @param observations         Sequence of observations
     * @param transitionMatrix     Provides the probabilities with which state transitions occur.
     *                             transitionMatrix[j][k] = 0.25 means that there is a chance of 0.25 to go from state
     *                             j to state k. Size must be |States| x |States| .
     * @param emissionMatrix       Provides the probabilities with which an emission occurs, given a state.
     *                             emissionMatrix[j][o] = 0.25 means that there is a 0.25 chance to emit an observation o
     *                             when in state j
     * @return the viterbi path
     */
    public static ViterbiResult calc(int[] observations, double[][] transitionMatrix, double[][] emissionMatrix) {
        int countStateSpace = transitionMatrix.length;
        double[][] viterbiVar = new double[countStateSpace][observations.length + 1];
        int[][] backtrackingVar = new int[countStateSpace][observations.length + 1];

        viterbiVar[0][0] = 1;
        
        /*
            Logarithm is applied element wise to the transition and emission matrices. Computing the logithms of
            the probabilities instead of the probabilities themselves ensures that the values stay in ranges which
            can be represented by doubles.
         */
        double[][] logTransitionMatrix = Viterbi.toLog(transitionMatrix.clone());
        double[][] logEmissionMatrix = Viterbi.toLog(emissionMatrix.clone());

        double maxProbability = calcViterbiBacktrackVars(observations, countStateSpace, viterbiVar, backtrackingVar, logTransitionMatrix, logEmissionMatrix);

        int[] path = reconstructOptimalPath(observations.length, countStateSpace, backtrackingVar, viterbiVar);

        return new ViterbiResult(path, maxProbability);
    }

    private static double calcViterbiBacktrackVars(int[] observations, int countStateSpace, double[][] viterbiVar, int[][] backtrackingVar, double[][] logTransitionMatrix, double[][] logEmissionMatrix) {
        double maxProbability = 0.;
        for (int observationIdx = 1; observationIdx <= observations.length; observationIdx++) {
            for (int state = 0; state < countStateSpace; state++) {
                double maxScore = viterbiVar[0][observationIdx - 1] + logTransitionMatrix[0][state];
                int argMaxScore = 0;
                for (int i = 1; i < countStateSpace; i++) {
                    double score = viterbiVar[i][observationIdx - 1] + logTransitionMatrix[i][state];
                    if (score > maxScore) {
                        maxScore = score;
                        argMaxScore = i;
                    }
                }
                maxScore += logEmissionMatrix[state][observations[observationIdx - 1]];
                viterbiVar[state][observationIdx] = maxScore;
                maxProbability = maxScore;
                backtrackingVar[state][observationIdx] = argMaxScore;
            }
        }
        return maxProbability;
    }

    private static int[] reconstructOptimalPath(int pathLength, int countStateSpace, int[][] backtrackingVar, double[][] viterbiVar) {
        int[] path = new int[pathLength];
        int last = path.length - 1;
        path[last] = 0;
        double maxLastColumn = viterbiVar[0][last + 1];
        for (int state = 1; state < countStateSpace; state++) {
            if (viterbiVar[state][last + 1] > maxLastColumn) {
                maxLastColumn = viterbiVar[state][last + 1];
                path[last] = state;
            }
        }
        for (int observation = last; observation > 0; observation--) {
            path[observation - 1] = backtrackingVar[path[observation]][observation + 1];
        }
        return path;
    }

    public static ViterbiResult calc(int[] observations, ProfileHMM profileHMM) {
        var transitionMatrix = profileHMM.getTransitionMatrix();
        var emissionMatrix = profileHMM.getEmissionMatrix();
        int countStateSpace = transitionMatrix.length;
        double[][] viterbiVar = new double[countStateSpace][observations.length + 2];
        int[][] backtrackingVar = new int[countStateSpace][observations.length + 2];

        viterbiVar[0][0] = 1;

        /*
            Logarithm is applied element wise to the transition and emission matrices. Computing the logithms of
            the probabilities instead of the probabilities themselves ensures that the values stay in ranges which
            can be represented by doubles.
         */
        double[][] logTransitionMatrix = Viterbi.toLog(transitionMatrix.clone());
        double[][] logEmissionMatrix = Viterbi.toLog(emissionMatrix.clone());

        double maxProbability = calcViterbiVars(observations, profileHMM, viterbiVar, backtrackingVar, logTransitionMatrix, logEmissionMatrix);
        var path = reconstructOptimalPath(observations.length, profileHMM, viterbiVar, backtrackingVar);
        return new ViterbiResult(path, maxProbability);
    }


    private static double calcViterbiVars(int[] observations, ProfileHMM profileHMM, double[][] viterbiVar, int[][] backtrackVars, double[][] logTransitionMatrix, double[][] logEmissionMatrix) {
        var countStateSpace = logTransitionMatrix.length;
        for (int observationIdx = 0; observationIdx < observations.length; observationIdx++) {
            // we begin at column 1;
            var viterbiIdx = observationIdx + 1;
            for (int state = 1; state < countStateSpace; state++) {
                var predecessorStates = profileHMM.getPossiblePredecessorIndeces(state);
                var compoundProbabilities = new ArrayList<Map.Entry<Integer, Double>>();
                for (int predecessor : predecessorStates) {
                    if (predecessor <= profileHMM.getLastInsert()) {
                        var entry = new AbstractMap.SimpleEntry<>(
                                predecessor,
                                logEmissionMatrix[state][observationIdx] +
                                        viterbiVar[predecessor][viterbiIdx - 1] +
                                        logTransitionMatrix[predecessor][state]
                        );
                        compoundProbabilities.add(entry);
                        entry = new AbstractMap.SimpleEntry<>(
                                predecessor,
                                logEmissionMatrix[state][observationIdx] +
                                        viterbiVar[predecessor][viterbiIdx - 1] +
                                        logTransitionMatrix[predecessor][state]
                        );
                        compoundProbabilities.add(entry);
                    } else {
                        var entry = new AbstractMap.SimpleEntry<>(
                                predecessor,
                                viterbiVar[predecessor][viterbiIdx] +
                                        logTransitionMatrix[predecessor][state]
                        );
                        compoundProbabilities.add(entry);
                    }
                }
                var argMaxAndMax = getMaxAndArgMax(compoundProbabilities);
                viterbiVar[state][viterbiIdx] = argMaxAndMax.getValue();
                backtrackVars[state][viterbiIdx] = argMaxAndMax.getKey();
            }
        }

        // termination
        var predecessorsToEndState = profileHMM.getPossiblePredecessorIndeces(profileHMM.getEndMatch());
        ArrayList<Map.Entry<Integer, Double>> lastColumnValsAndPredStates = predecessorsToEndState
                .stream()
                .map(state -> new AbstractMap.SimpleEntry<>(state, viterbiVar[state][viterbiVar.length-1] + logTransitionMatrix[state][profileHMM.getEndMatch()]))
                .collect(Collectors.toCollection(ArrayList::new));
        var argMaxAndMax = getMaxAndArgMax(lastColumnValsAndPredStates);
        viterbiVar[profileHMM.getEndMatch()][viterbiVar.length - 1] = argMaxAndMax.getValue();
        backtrackVars[profileHMM.getEndMatch()][backtrackVars.length - 1] = argMaxAndMax.getKey();
        return argMaxAndMax.getValue();
    }

    private static int[] reconstructOptimalPath(int observationCount, ProfileHMM profileHMM, double[][] viterbiVars, int[][] backtrackingVar) {
        var path = new int[observationCount + 2]; // for end and begin state
        path[path.length-1] = profileHMM.getEndMatch();

        for (int observation = observationCount-1; observation >= -1; observation--) {
            path[observation+1] = backtrackingVar[path[observation+2]][observation + 2];
        }
        return path;
    }

    /**
     * Transforms the provided matrix by an element wise natural logarithm in place.
     *
     * @param matrix 2d matrix that should be transformed.
     * @return the transformed matrix.
     */
    private static double[][] toLog(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = Math.log(matrix[i][j]);
            }
        }
        return matrix;
    }

    private static Map.Entry<Integer, Double> getMaxAndArgMax(List<Map.Entry<Integer, Double>> input) {
        if (input.isEmpty()) {
            throw new IllegalArgumentException("Can not call getMaxAndArgMax with empty List");
        }
        var iter = input.iterator();
        var max = iter.next();

        while (iter.hasNext()) {
            var next = iter.next();
            if (next.getValue() >= max.getValue()) {
                max = next;
            }
        }
        return max;
    }
}
