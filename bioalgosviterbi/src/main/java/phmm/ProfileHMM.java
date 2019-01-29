package phmm;

import fasta.FASTASequence;
import util.Util;

import java.util.*;

public class ProfileHMM {
    private double[][] transitionMatrix;
    private double[][] emissionMatrix;

    final private char gapSymbol;
    final private double pseudoCount;
    final private int observationStatesCount;
    final private int columnCount;

    private int beginMatch;
    private int endMatch;
    private int firstInsert;
    private int lastInsert;
    private int firstDelete;
    private int lastDelete;

    private int stateCount;

    private boolean appliedPCAndNormalized = false;

    public ProfileHMM(List<FASTASequence> sequences,
                      char gapSymbol,
                      Map<Character, Integer> observationMap,
                      double pseudocount,
                      double matchThreshold) {
        if (sequences.size() == 0 || observationMap.size() == 0) {
            throw new IllegalArgumentException("Input sequences and observationMap cannot be empty");
        } else if (pseudocount < 0) {
            throw new IllegalArgumentException("No negative pseudocounts permitted");
        } else if (observationMap.containsKey(gapSymbol)) {
            throw new IllegalArgumentException("Obersvation map cannot contain gap symbol");
        }

        this.gapSymbol = gapSymbol;
        this.pseudoCount = pseudocount;
        this.observationStatesCount = observationMap.size();
        this.columnCount = sequences.get(0).getSequence().length;

        var matchColumns = getMatchColumns(sequences, matchThreshold);

        calcStateBorders(matchColumns.size());
        calcEmissionMatrix(sequences, observationMap, matchColumns.clone());
        calcTransitionMatrix(sequences, matchColumns);

        Util.toLog(transitionMatrix);
        Util.toLog(emissionMatrix);
    }

    public double[][] getTransitionMatrix() {
        return transitionMatrix;
    }

    public double[][] getEmissionMatrix() {
        return emissionMatrix;
    }

    public int getBeginMatch() {
        return beginMatch;
    }

    public int getEndMatch() {
        return endMatch;
    }

    public int getFirstInsert() {
        return firstInsert;
    }

    public int getLastInsert() {
        return lastInsert;
    }

    public int getFirstDelete() {
        return firstDelete;
    }

    public int getLastDelete() {
        return lastDelete;
    }

    public int getStateCount() {
        return stateCount;
    }

    public int stateToColum(int state) {
        if (state < 0 || state > lastDelete) {
            throw new IllegalArgumentException("Can not call stateToColum with < 0 or > lastDelete");
        }
        if (state <= endMatch) {
            return state;
        } else if (state <= lastInsert) {
            return state - firstInsert;
        } else {
            return state - firstDelete +1;
        }
    }

    public int nextKindOfState(int state) {
        if (state < 0 || state > lastDelete) {
            throw new IllegalArgumentException("Can not call nextKindOfState with < 0 or > lastDelete");
        }
        if (state <= endMatch) {
            return firstInsert;
        } else if (state <= lastInsert) {
            return firstDelete;
        } else {
            return lastDelete;
        }
    }

    public List<Integer> getPossibleSuccessorIndeces(int index) {
        var successors = new ArrayList<Integer>();
        if (index < 0 || index > lastDelete) {
            throw new IllegalArgumentException("Invalid index, must be >= 0, <= lastDelete State");
        }

        var insertCount = lastInsert - firstInsert + 1;

        if (index < endMatch) { // skipping end match state
            successors.add(index+1); // match state
            successors.add(firstInsert+index); // insert state
            // last match and end match state have no transition to delete
            if (index < endMatch -1) {
                successors.add(firstDelete+index);
            }
        } else if (index >= firstInsert && index <= lastInsert) {
            successors.add(index - firstInsert + 1); // match state
            successors.add(index);  // insert
            if (index < lastInsert) {
                successors.add(insertCount + index); // delete
            }
        } else if (index >= firstDelete){
            successors.add(index - lastInsert +1); // match
            successors.add(index - insertCount); // insert
            if (index < lastDelete) {
                successors.add(index+1);    // delete
            }
        }

        return successors;
    }

    public List<Integer> getPossiblePredecessorIndeces(int index) {
        if (index < 0 || index > lastDelete) {
            throw new IllegalArgumentException("No predecessor states for out of range index");
        }
        var predecessors = new ArrayList<Integer>();
        if (index > beginMatch && index <= endMatch) {
            predecessors.add(index-1);      // prev match
            predecessors.add(firstInsert + index -1);   // prev insert
            if (index > beginMatch+1) {
                predecessors.add(firstDelete + index -2);
            }
        } else if (index <= lastInsert) {
            predecessors.add(index - firstInsert);     //prev match
            predecessors.add(index);        // prev inser
            if (index > firstInsert) {
                predecessors.add(index - firstInsert + firstDelete -1);
            }
        } else {
            predecessors.add(index - firstDelete);  // prev match
            predecessors.add(index - firstDelete + firstInsert);    // prev insert
            if (index > firstDelete) {
                predecessors.add(index - 1);
            }
        }

        return predecessors;
    }

    private ArrayDeque<Integer> getMatchColumns(List<FASTASequence> sequences, double matchThreshold) {
        var matchColums = new ArrayDeque<Integer>();
        var columnCount = sequences.get(0).getSequence().length;

        for (var column=0; column<columnCount; column++) {
            var gapCount = 0;
            for (FASTASequence seq : sequences) {
                if (seq.getSequence()[column] == gapSymbol) {
                    gapCount++;
                }
            }
            if (gapCount < sequences.size() * matchThreshold) {
                matchColums.add(column);
            }
        }
        return matchColums;
    }

    private void calcStateBorders(int matchStateWithoutBeginningAndEnd) {
        var matchCount = matchStateWithoutBeginningAndEnd;
        beginMatch = 0;
        endMatch = matchCount + 1;
        firstInsert = endMatch + 1;
        lastInsert = firstInsert + matchCount;
        firstDelete = lastInsert + 1;
        lastDelete = firstDelete + matchCount - 1;
        stateCount = lastDelete +1;
    }

   private void calcEmissionMatrix(List<FASTASequence> sequences,
                                   Map<Character, Integer> observationMap,
                                   Deque<Integer> matchColumns) {
        // deletes have no emissions
        emissionMatrix = new double[lastInsert + 1][observationStatesCount];

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
                    emissionMatrix[firstInsert + profileColumn][mappedElem]++;
                }
            }
        }

        // skip begin state row
        for (var i=1; i<emissionMatrix.length; i++) {
            if (i == endMatch)
                continue;       // end match state has no emissions
            var rowsum = Arrays.stream(emissionMatrix[i]).sum();
            var divisor = rowsum + pseudoCount * observationStatesCount;
            for(int j=0; j<emissionMatrix[i].length; j++) {
                emissionMatrix[i][j] = (emissionMatrix[i][j] + pseudoCount) / divisor;
            }
        }
    }

    private void calcTransitionMatrix(List<FASTASequence> sequences, Deque<Integer> matchColumns) {
        transitionMatrix = new double[stateCount][stateCount];

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
            for (FASTASequence FASTASeq : sequences) {
                var seq = FASTASeq.getSequence();
                if (!isMatchColumn && seq[column] == gapSymbol) {
                    continue;
                }

                var prevState = getPreviousState(column,
                        seq,
                        isMatchColumn ? previousMatchColumn : currentMatchColumn,
                        isMatchColumn ? profileColumn - 1 : profileColumn);

                var currState = profileColumn;

                if (!isMatchColumn) {   // we are in an insert state
                    currState += firstInsert;
                } else if (seq[column] == gapSymbol) {      // we are in a delete state
                     currState += firstDelete - 1;      // deletes start in column 1
                }

                transitionMatrix[prevState][currState]++;
            }
        }

        // handle transitions to end state
        var lastColumn = columnCount -1;
        for (var FASTASeq : sequences) {
            var seq = FASTASeq.getSequence();
            var prevState = getPreviousState(lastColumn + 1,
                    seq,
                    currentMatchColumn ,
                    profileColumn);
            transitionMatrix[prevState][endMatch]++;
        }

        applyPseudcountAndNormalizeToTransitionMatrix();
    }

    private int getPreviousState(int index,
                                 char[] sequence,
                                 Optional<Integer> previousMatchColumnFromIndex,
                                 int prevProfileColumn) {
        var previousMatchColumn = previousMatchColumnFromIndex.orElse(-1);

        for (var i=index-1; i>=previousMatchColumnFromIndex.orElse(-1); i--) {
            // if we do not find a previous state from the index given the sequence
            // then the only one left is the beginning state.
            // we have to look at all the earlier observations first, because the previous state
            // could also be the first insert
            if (i == -1) {
                return beginMatch;
            } else if (i == previousMatchColumn) {
                if (sequence[i] == gapSymbol) {
                    //first delete is in column 1
                    return firstDelete + prevProfileColumn -1;
                } else {
                    return beginMatch + prevProfileColumn;
                }
            } else {
                if (sequence[i] == gapSymbol)
                    continue;
                return firstInsert + prevProfileColumn;
            }
        }
        throw new IllegalArgumentException("Error in previous state.");
    }

    private void applyPseudcountAndNormalizeToTransitionMatrix() {
        if (this.appliedPCAndNormalized) {
            throw new IllegalStateException("Can not call applyPseudcountAndNormalizeToTransitionMatrix twice!");
        }
        for(var fromState=0; fromState < transitionMatrix.length; fromState++) {
            for (int toState : getPossibleSuccessorIndeces(fromState)) {
                transitionMatrix[fromState][toState] += pseudoCount;
            }
        }
        for(var fromState=0; fromState < transitionMatrix.length; fromState++) {
            if (fromState == endMatch) // has no transitions
                continue;
            var rowSum = Arrays.stream(transitionMatrix[fromState]).sum();
            for (var toState=0; toState<transitionMatrix.length; toState++) {
                transitionMatrix[fromState][toState] /= rowSum;
            }
        }
    }


}

