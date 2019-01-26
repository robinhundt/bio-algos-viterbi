package viterbi;

import static org.junit.Assert.assertArrayEquals;

import fasta.FASTASequence;
import org.junit.Test;
import phmm.ProfileHMM;
import util.Util;

import java.util.ArrayList;

public class ViterbiTest {
    @Test
    public void testViterbiWithCasinoData() {
        int[] observations = {2, 0, 4, 0, 0, 5, 1, 3, 5, 3, 3, 5, 5, 3, 3, 1, 3, 4, 2, 0, 0, 2, 1, 0, 5, 2, 0, 0, 5, 3, 0, 4, 1, 0, 2, 2, 5, 1, 4, 0, 3, 3, 4, 3, 2, 5, 2, 0, 5, 4, 5, 5, 1, 5, 4, 5, 5, 5, 5, 5, 5, 4, 0, 0, 5, 5, 3, 4, 2, 0, 2, 1, 5, 4, 0, 1, 3, 4, 5, 2, 5, 5, 5, 3, 5, 2, 0, 5, 2, 5, 5, 5, 2, 0, 5, 1, 2, 1, 5, 3, 4, 4, 1, 2, 5, 1, 5, 5, 5, 5, 5, 5, 1, 4, 0, 4, 0, 5, 2, 0, 1, 1, 1, 4, 4, 4, 3, 3, 0, 5, 5, 5, 4, 5, 5, 4, 5, 2, 4, 5, 3, 2, 1, 3, 2, 5, 3, 0, 2, 0, 4, 0, 2, 3, 5, 4, 0, 3, 5, 2, 4, 2, 3, 0, 0, 0, 1, 5, 3, 0, 3, 5, 1, 5, 1, 4, 2, 2, 4, 5, 2, 5, 5, 0, 5, 2, 5, 5, 5, 3, 5, 5, 1, 2, 1, 4, 2, 3, 3, 0, 2, 5, 5, 0, 5, 5, 0, 0, 5, 2, 1, 4, 1, 4, 5, 1, 3, 5, 1, 1, 4, 4, 1, 5, 4, 1, 4, 1, 1, 5, 5, 3, 2, 4, 2, 4, 2, 2, 2, 5, 1, 2, 2, 0, 1, 0, 5, 1, 4, 2, 5, 3, 3, 0, 3, 3, 2, 1, 2, 2, 4, 0, 5, 2, 1, 3, 2, 5, 2, 2, 5, 5, 4, 4, 5, 1, 3, 5, 5, 5, 5, 1, 5, 2, 1, 5, 5, 5, 5, 0, 1, 2, 4, 4, 1, 3, 4, 1, 3, 1};

        /*
             S      F      L
            S 0     0.5    0.5
            F 0     0.95  0.05
            L 0     0.1   0.9
         */
        double[][] transitionMatrix = {
                {0, 0.5, 0.5},
                {0, 0.95, 0.05},
                {0, 0.1, 0.9}
        };
        /*
             1   2   3   ...
           S 
           F probabilities
           L
         */
        double[][] emissionMatrix = {
                {0., 0., 0., 0., 0., 0.},
                {1./6, 1./6, 1./6, 1./6, 1./6, 1./6},
                {1./10, 1./10, 1./10, 1./10, 1./10, 1./2}
        };

        int[] expectedOutput = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

        ViterbiResult viterbiResult = Viterbi.calc(observations, transitionMatrix, emissionMatrix);
        System.out.println(viterbiResult.getMaxProbability());
        assertArrayEquals(expectedOutput, viterbiResult.getViterbiPath());
    }

    @Test
    public void testViterbiPHMM() {
        final var pseudoCount = 0.00001;

        var sequences = new ArrayList<FASTASequence>();

        sequences.add(new FASTASequence("1", new char[]{'-', 'T', 'G'}));
        sequences.add(new FASTASequence("2", new char[]{'A', '-', 'G'}));
        sequences.add(new FASTASequence("3", new char[]{'A', '-', 'G'}));

        var observationMap = Util.createObersavtionMap();

        var profileHmm = new ProfileHMM(sequences, '-', observationMap, pseudoCount, 0.5);

        var observations = new int[]{3, 3, 3, 2};

        var viterbiResult = Viterbi.calc(observations, profileHmm);


    }
}
