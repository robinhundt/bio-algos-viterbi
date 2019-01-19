package phmm;

import fasta.FASTASequence;
import static org.junit.Assert.assertArrayEquals;
import org.junit.Test;


import java.util.ArrayList;
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

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 0);

        var expectedEmissionMatrixCount = new double[][]{
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

        for (var i=1; i<expectedEmissionMatrixCount.length; i++) {
            if (i==5)
                continue;
            for (var j=; j<expectedEmissionMatrixCount[0].length) {

            }
        }



        var expectedEmissionMatrix = expectedEmissionMatrixCount.clone();

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

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 0);

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

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, 0);

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

        assertArrayEquals(expectedEmissionMatrix, profileHmm.getEmissionMatrix());

    }

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
