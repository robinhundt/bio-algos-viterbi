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

        calcEmissionMatrix(sequences, gapSymbol, observationMap, pseudoCount, matchColumns.clone(), columnCount, observationStatesCount, matchStatesCount, insertionStatesCount, deletionStatesCount, statesCount);
        calcTransitionMatrix(sequences, gapSymbol, observationMap, pseudoCount, matchColumns, columnCount, observationStatesCount, matchStatesCount, insertionStatesCount, deletionStatesCount, statesCount);
    }


    private void calcTransitionMatrix(ArrayList<FASTASequence> sequences, Character gapSymbol, HashMap<Character, Integer> observationMap, int pseudoCount, ArrayDeque<Integer> matchColumns, int columnCount, int observationStatesCount, int matchStatesCount, int insertionStatesCount, int deletionStatesCount, int statesCount) {
        transitionMatrix = new double[statesCount][statesCount];

        Optional<Integer> previousMatchColumn = Optional.empty();
        Optional<Integer> currentMatchColumn = Optional.empty();
        var profileColumn = 0;


        for (var column=0; column<columnCount; column++) {
            if (!matchColumns.isEmpty() && column == matchColumns.peek()) {
                previousMatchColumn = currentMatchColumn;
                currentMatchColumn = Optional.of(matchColumns.removeFirst());
                profileColumn++;
            }

            var isMatchColumn = currentMatchColumn.equals(Optional.of(column));

            for (FASTASequence FASTAseq : sequences) {
                var seq = FASTAseq.getSequence();

                if (column == 0) {
                    if (isMatchColumn && seq[column] == gapSymbol) {
                        transitionMatrix[0][matchStatesCount+insertionStatesCount]++;
                    } else {
                        if (isMatchColumn) {
                            transitionMatrix[0][1]++;
                        } else if (seq[column] != gapSymbol) {
                            transitionMatrix[0][matchStatesCount]++;
                        }
                    }
                } else if (isMatchColumn || seq[column] != gapSymbol){
                    var prevState = getPreviousState(profileColumn,
                            seq,
                            column,
                            isMatchColumn ? profileColumn - 1 : profileColumn,
                            isMatchColumn ? previousMatchColumn.orElse(-1) : currentMatchColumn.orElse(-1),
                            gapSymbol, isMatchColumn, matchStatesCount, insertionStatesCount);

                    var currState = profileColumn;
                    if (isMatchColumn && seq[column] == gapSymbol) {
                        currState += matchStatesCount + insertionStatesCount -1;
                    } else if (!isMatchColumn) {
                        currState += matchStatesCount;
                    }

                    transitionMatrix[prevState][currState]++;
                }


            }
        }

        var lastColumn = columnCount-1;
        var isMatchColumn = currentMatchColumn.equals(Optional.of(lastColumn));
        // handle transtitions to end state
        for (var FASTAseq : sequences) {
            var seq = FASTAseq.getSequence();
            if (isMatchColumn) {
                var fromState = seq[lastColumn] == gapSymbol ? statesCount -1 : matchStatesCount -2;
                transitionMatrix[fromState][matchStatesCount-1]++;
            } else {
                int prevState = matchStatesCount+insertionStatesCount-1;
                if (seq[lastColumn] == gapSymbol) {
                    prevState = getPreviousState(profileColumn,
                            seq,
                            lastColumn,
                            profileColumn,
                            currentMatchColumn.orElse(-1),
                            gapSymbol, isMatchColumn, matchStatesCount, insertionStatesCount);
                }
                        transitionMatrix[prevState][matchStatesCount-1]++;
            }
        }

        applyPseudcountAndNormalize(transitionMatrix, pseudoCount, matchStatesCount, insertionStatesCount, deletionStatesCount);
    }

    private static void applyPseudcountAndNormalize(double[][] transitionMatrix, int pseudocount, int matchCount, int insertCount, int deleteCount) {
        for(var fromState=0; fromState < transitionMatrix.length; fromState++) {
            for (int toState : getPossibleSuccessorIndeces(fromState, matchCount, insertCount, deleteCount)) {
                transitionMatrix[fromState][toState] += pseudocount;
            }
        }
        for(var fromState=0; fromState < transitionMatrix.length; fromState++) {
            var rowSum = Arrays.stream(transitionMatrix[fromState]).sum();
            if (rowSum == 0)
                continue;
            for (var toState=0; toState<transitionMatrix.length; toState++) {
                transitionMatrix[fromState][toState] /= rowSum;
            }
        }
    }

    private int getPreviousState(int currProfileColumn, char[] sequence, int index,
                                 int previousMatchState, int previousMatchColumn, char gapSymbol,
                                 boolean fromMatchState, int matchStatesCount, int insertStatesCount) {

        for (var i=index-1; i>=previousMatchColumn; i--) {
            if (i == -1) {
                return 0;
            } else if (i == previousMatchColumn) {
                if (sequence[i] == gapSymbol) {
                    var deleteState = currProfileColumn + matchStatesCount + insertStatesCount - 1;
                    return fromMatchState ? deleteState -1 : deleteState;
                } else {
                    return previousMatchState;
                }
            } else {
                if (sequence[i] == gapSymbol)
                    continue;
                var insertState = currProfileColumn + matchStatesCount;
                return fromMatchState ? insertState - 1 : insertState;
            }
        }
        throw new IllegalArgumentException("Error in previous state.");
    }


    private void calcEmissionMatrix(ArrayList<FASTASequence> sequences, Character gapSymbol, HashMap<Character, Integer> observationMap, int pseudoCount, ArrayDeque<Integer> matchColumns, int columnCount, int observationStatesCount, int matchStatesCount, int insertionStatesCount, int deletionStatesCount, int statesCount) {
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

    public static ArrayList<Integer> getPossibleSuccessorIndeces(int index, int matchCount, int insertCount, int deleteCount) {
        var successors = new ArrayList<Integer>();
        if (index >= 0 && index < matchCount-1) { // skipping end match state
            successors.add(index+1); // match state
            successors.add(matchCount+index); // insert state
            // last match and end match state have no transition to delete
            if (index < matchCount-2) {
                successors.add(matchCount+insertCount+index);
            }
        } else if (index >= matchCount && index < matchCount + insertCount) {
            successors.add(index - matchCount + 1); // match state
            successors.add(index);  // insert
            if (index < matchCount + insertCount - 1) {
                successors.add(insertCount + index); // delete
            }
        } else if (index >= matchCount + insertCount && index < matchCount+insertCount+deleteCount){
            successors.add(index - insertCount);
            successors.add(index - insertCount - matchCount +2);
            if (index < matchCount + insertCount + deleteCount -1) {
                successors.add(index+1);
            }
        }

        return successors;
    }

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
