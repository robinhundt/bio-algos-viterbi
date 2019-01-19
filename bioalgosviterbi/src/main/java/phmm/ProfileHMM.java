package phmm;

import fasta.FASTASequence;

import java.util.*;

public class ProfileHMM {

    private double[][] transitionMatrix;
    private double [][] emissionMatrix;

    public ProfileHMM(ArrayList<FASTASequence> sequences,
                      Character gapSymbol,
                      HashMap<Character, Integer> observationMap,
                      int pseudoCount) {
        var matchColumns = getMatchColumns(sequences, gapSymbol);
        var columnCount = sequences.get(0).getSequence().length;

        if (matchColumns.size() == 0) {
            throw new IllegalArgumentException("No match columns!");
        }

        final var observationStatesCount = observationMap.size();
        final var matchStatesCount = matchColumns.size() + 2; // add begin and end
        final var insertionStatesCount = matchColumns.size() + 1;
        final var deletionStatesCount = matchColumns.size();
        final var statesCount = matchStatesCount + insertionStatesCount + deletionStatesCount;

        transitionMatrix = new double[statesCount][statesCount];
        emissionMatrix = new double[statesCount][observationStatesCount];

        Optional<Integer> currentMatchColumn = Optional.empty();
        var profileColumn = 0;


        for (var column=0; column<columnCount; column++) {
            if (!matchColumns.isEmpty() && column == matchColumns.peek()) {
                currentMatchColumn = Optional.of(matchColumns.removeFirst());
                profileColumn++;
            }

            var isMatchColumn = currentMatchColumn.equals(Optional.of(column));

            for (FASTASequence seq : sequences) {
                var elem = seq.getSequence()[column];
                if (elem == gapSymbol)
                    continue;

                var mappedElem = observationMap.get(elem);

                if (isMatchColumn) {
                    emissionMatrix[profileColumn][mappedElem]++;
                } else {
                    emissionMatrix[matchStatesCount + profileColumn][mappedElem]++;
                }
            }
        }

        // skip begin state row
        for (var i=1; i<matchStatesCount+insertionStatesCount; i++) {
            var rowsum = Arrays.stream(emissionMatrix[i]).sum();
            var divisor = rowsum + pseudoCount * (observationStatesCount - 1);
            for(int j=0; j<emissionMatrix[i].length-1; j++) {
                emissionMatrix[i][j] = (emissionMatrix[i][j] + pseudoCount) / divisor;
            }
        }
        // end state has no emissions
        for (var i=0; i < emissionMatrix[matchStatesCount-1].length; i++) {
            emissionMatrix[matchStatesCount-1][i] = 0;
        }

        for (var row=statesCount-deletionStatesCount; row<statesCount; row++) {
            emissionMatrix[row][observationMap.get(gapSymbol)] = 1;
        }


    }

    public double[][] getTransitionMatrix() {
        return transitionMatrix;
    }

    public double[][] getEmissionMatrix() {
        return emissionMatrix;
    }

//    public ArrayList<int[]> getPossibleSuccessorIndeces(int index, int matchCount, int insertCount, int deleteCount) {
//        var successors = new ArrayList<int[]>();
//        if (index < matchCount) {
//
//        }
//    }

    private ArrayDeque<Integer> getMatchColumns(ArrayList<FASTASequence> sequences, Character gapSymbol) {
        var matchColums = new ArrayDeque<Integer>();
        var columnCount = sequences.get(0).getSequence().length;

        for (var column=0; column<columnCount; column++) {
            var gapCount = 0;
            for (FASTASequence seq : sequences) {
                if (seq.getSequence()[column] == gapSymbol) {
                    gapCount++;
                }
            }
            if (gapCount < sequences.size() / 2.0) {
                matchColums.add(column);
            }
        }
        return matchColums;
    }

    private static Integer getMappedNthSequenceElem(FASTASequence seq,int n, HashMap<Character, Integer> observationMap) {
        return observationMap.get(seq.getSequence()[n]);
    }
}
