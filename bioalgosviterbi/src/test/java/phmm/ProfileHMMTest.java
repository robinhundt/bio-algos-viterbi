package phmm;

import fasta.FASTASequence;
import static org.junit.Assert.assertArrayEquals;
import org.junit.Test;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class ProfileHMMTest {

    @Test
    public void testEmissionMatrixGenerationFirstColumnIsMatch() {
        var sequences = new ArrayList<FASTASequence>();

        sequences.add(new FASTASequence("1", new char[]{'A', 'G', '-', '-', '-', 'C'}));
        sequences.add(new FASTASequence("2", new char[]{'A', '-', 'A', 'G', '-', 'C'}));
        sequences.add(new FASTASequence("3", new char[]{'A', 'G', '-', 'A', 'A', '-'}));
        sequences.add(new FASTASequence("4", new char[]{'-', '-', 'A', 'A', 'A', 'C'}));
        sequences.add(new FASTASequence("5", new char[]{'A', 'G', '-', '-', '-', 'C'}));

        var observationMap = createObersavtionMap();

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 1);

        var expectedEmissionMatrix = new double[][]{
                {0, 0, 0, 0, 0},   // M0
                {4, 0, 0, 0, 0},
                {0, 0, 3, 0, 0},
                {2, 0, 1, 0, 0},
                {0, 4, 0, 0, 0},
                {0, 0, 0, 0, 0},   // M5
                {0, 0, 0, 0, 0},   // I0
                {0, 0, 0, 0, 0},
                {2, 0, 0, 0, 0},
                {2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0},   // I4
                {0, 0, 0, 0, 1},   // D0
                {0, 0, 0, 0, 1},
                {0, 0, 0, 0, 1},
                {0, 0, 0, 0, 1},   // D3
        };

        //skip first and last match state, also delete states
        for (var i=1; i<expectedEmissionMatrix.length - 4; i++) {
            if (i==5)
                continue;
            var rowsum = Arrays.stream(expectedEmissionMatrix[i]).sum();
            var divisor = rowsum + 4; // add pseudocount
            for(int j=0; j<expectedEmissionMatrix[i].length -1; j++) {
                expectedEmissionMatrix[i][j] = (expectedEmissionMatrix[i][j] + 1) / divisor;
            }
        }

        assertArrayEquals(expectedEmissionMatrix, profileHmm.getEmissionMatrix());

    }

    @Test
    public void testEmissionMatrixGenerationFirstColumnIsNoMatch() {
        var sequences = new ArrayList<FASTASequence>();

        sequences.add(new FASTASequence("1", new char[]{'-', 'G', '-', '-', '-', 'C'}));
        sequences.add(new FASTASequence("2", new char[]{'-', '-', 'A', 'G', '-', 'C'}));
        sequences.add(new FASTASequence("3", new char[]{'A', 'G', '-', 'A', 'A', '-'}));
        sequences.add(new FASTASequence("4", new char[]{'-', '-', 'A', 'A', 'A', 'C'}));
        sequences.add(new FASTASequence("5", new char[]{'A', 'G', '-', '-', '-', 'C'}));

        var observationMap = createObersavtionMap();

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 1);

        var expectedEmissionMatrix = new double[][]{
                {0, 0, 0, 0, 0},   // M0
                {0, 0, 3, 0, 0},
                {2, 0, 1, 0, 0},
                {0, 4, 0, 0, 0},
                {0, 0, 0, 0, 0},   // M4
                {2, 0, 0, 0, 0},   // I0
                {2, 0, 0, 0, 0},
                {2, 0, 0, 0, 0},
                {0, 0, 0, 0, 0},   // I3
                {0, 0, 0, 0, 1},   // D0
                {0, 0, 0, 0, 1},
                {0, 0, 0, 0, 1},   // D2
        };

        //skip first and last match state, also delete states
        for (var i=1; i<expectedEmissionMatrix.length - 3; i++) {
            if (i==4)
                continue;
            var rowsum = Arrays.stream(expectedEmissionMatrix[i]).sum();
            var divisor = rowsum + 4; // add pseudocount
            for(int j=0; j<expectedEmissionMatrix[i].length -1; j++) {
                expectedEmissionMatrix[i][j] = (expectedEmissionMatrix[i][j] + 1) / divisor;
            }
        }

        assertArrayEquals(expectedEmissionMatrix, profileHmm.getEmissionMatrix());

    }

    @Test
    public void testEmissionMatrixGenerationLastColumnIsNoMatch() {
        var sequences = new ArrayList<FASTASequence>();

        sequences.add(new FASTASequence("1", new char[]{'-', 'G', '-', '-', '-', '-'}));
        sequences.add(new FASTASequence("2", new char[]{'-', '-', 'A', 'G', '-', '-'}));
        sequences.add(new FASTASequence("3", new char[]{'A', 'G', '-', 'A', 'A', '-'}));
        sequences.add(new FASTASequence("4", new char[]{'-', '-', 'A', 'A', 'A', 'C'}));
        sequences.add(new FASTASequence("5", new char[]{'A', 'G', '-', '-', '-', 'C'}));

        var observationMap = createObersavtionMap();

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 1);

        var expectedEmissionMatrix = new double[][]{
                {0, 0, 0, 0, 0},   // M0
                {0, 0, 3, 0, 0},
                {2, 0, 1, 0, 0},
                {0, 0, 0, 0, 0},   // M4
                {2, 0, 0, 0, 0},   // I0
                {2, 0, 0, 0, 0},
                {2, 2, 0, 0, 0},   // I4
                {0, 0, 0, 0, 1},   // D0
                {0, 0, 0, 0, 1},   // D1
        };

        //skip first and last match state, also delete states
        for (var i=1; i<expectedEmissionMatrix.length - 2; i++) {
            if (i==3)
                continue;
            var rowsum = Arrays.stream(expectedEmissionMatrix[i]).sum();
            var divisor = rowsum + 4; // add pseudocount
            for(int j=0; j<expectedEmissionMatrix[i].length -1; j++) {
                expectedEmissionMatrix[i][j] = (expectedEmissionMatrix[i][j] + 1) / divisor;
            }
        }

        assertArrayEquals(expectedEmissionMatrix, profileHmm.getEmissionMatrix());

    }

//    @Test
//    public void testTransitionMatrixOnlyMatch() {
//        var sequences = new ArrayList<FASTASequence>();
//
//        sequences.add(new FASTASequence("1", new char[]{'A', '-', 'C'}));
//        sequences.add(new FASTASequence("2", new char[]{'A', 'G', 'C'}));
//        sequences.add(new FASTASequence("3", new char[]{'A', 'A', '-'}));
//        sequences.add(new FASTASequence("4", new char[]{'-', 'A', 'C'}));
//        sequences.add(new FASTASequence("5", new char[]{'A', '-', 'C'}));
//
//        var observationMap = createObersavtionMap();
//
//        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 1);
//
//        var transitionMatrix = new double[][]{
////                    0 1 2 3 4 0 1 2 3 1 2 3
///*0*/                {0,1,0,0,0,0,0,0,0,0,0,0},
///*1*/                {0,0,2,0,0,0,0,0,0,0,2,0},
///*2*/                {0,0,0,2,0,0,0,0,0,0,0,1},
///*3*/                {0,0,0,0,4,0,0,0,0,0,0,0},
///*4*/                {0,0,0,0,0,0,0,0,0,0,0,0},
///*0*/                {0,0,0,0,0,0,0,0,0,0,0,0},
///*1*/                {0,0,0,0,0,0,0,0,0,0,0,0},
///*2*/                {0,0,0,0,0,0,0,0,0,0,0,0},
///*3*/                {0,0,0,0,0,0,0,0,0,0,0,0},
///*1*/                {0,0,1,0,0,0,0,0,0,0,0,0},
///*2*/                {0,0,0,2,0,0,0,0,0,0,0,0},
///*3*/                {0,0,0,0,0,0,0,0,0,0,0,0},
//        };
//
//        //skip first and last match state, also delete states
//        for (var i=1; i<expectedEmissionMatrix.length - 4; i++) {
//            if (i==5)
//                continue;
//            var rowsum = Arrays.stream(expectedEmissionMatrix[i]).sum();
//            var divisor = rowsum + 4; // add pseudocount
//            for(int j=0; j<expectedEmissionMatrix[i].length -1; j++) {
//                expectedEmissionMatrix[i][j] = (expectedEmissionMatrix[i][j] + 1) / divisor;
//            }
//        }
//
//    }




    private HashMap<Character, Integer> createObersavtionMap() {
        var observationMap = new HashMap<Character, Integer>();
        observationMap.put('A', 0);
        observationMap.put('C', 1);
        observationMap.put('G', 2);
        observationMap.put('T', 3);
        observationMap.put('-', 4);
        return observationMap;
    }

}
