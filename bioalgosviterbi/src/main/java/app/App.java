package app;

import fasta.FASTAParser;
import phmm.ProfileHMM;
import viterbi.EmissionMatrix;
import viterbi.TransitionMatrix;
import viterbi.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.io.File;
import java.util.*;

/**
 * Hello world!
 *
 */
public class App 
{

    private static ArrayList<String> getFileList(String path) {
        ArrayList<String> testFiles = new ArrayList<String>();
        Path filePath = Paths.get(path);

        if (Files.isDirectory(filePath)) {
            File directory = new File(path);

            for (var file : directory.listFiles()) {
                if (file.isFile()) {
                    testFiles.add(file.getAbsolutePath());
                } else if (file.isDirectory()) {
                    testFiles.addAll(getFileList(file.getAbsolutePath()));
                }
            }
        } else if (Files.isRegularFile(filePath)) {
            testFiles.add(path);
        }

        return testFiles;
    }

    public static void main( String[] args )
    {
        if (args.length < 1) {
            System.out.println("Helpmessage");
            return;
        }
        try {
            Character gapSymbol = '-';
            HashMap<Character, Integer> observationMap = new HashMap<Character, Integer>(Map.of('A', 0, 'C', 1, 'G', 2, 'T', 3, 'U', 3, gapSymbol, 4));
            int pseudoCount = 1;
            // import training sequences
            var trainingSequences = FASTAParser.parse(Paths.get(args[0]));

            // create profil HMM with test sequences
            ProfileHMM pHMM = new ProfileHMM(trainingSequences, gapSymbol, observationMap, pseudoCount);
            double[][] emissionMatrix = pHMM.getEmissionMatrix();
            double[][] transitionMatrix = pHMM.getTransitionMatrix();


            //get test sequences
            ArrayList<String> testFiles = getFileList(args[1]);

            for (var testFile : testFiles) {
                var sequences = FASTAParser.parse(Paths.get(testFile));
                ArrayList<int[]> observationLines = new ArrayList<int[]>(sequences.size());
                ArrayList<String> vitProbabilities = new ArrayList<String>(testFile.length());

                for (var sequence : sequences) {
                    observationLines.add(sequence.parseBasesToInt(observationMap));
                }

                for (var observations: observationLines) {
                    //calculate viterbi path and probability
                    
                    ViterbiResult viterbiResult = Viterbi.calc(observations, transitionMatrix, emissionMatrix);
                    vitProbabilities.add(String.valueOf(viterbiResult.getMaxProbability()));
                    
                    ;
                }
                // create new file and store each probability in one line
                /*
                File results = new File(testFile);
                results.createNewFile();
                Files.write(results.getAbsolutePath(), vitProbabilities);;
                */

                // System.out.println(observationLines.get(0));
                System.out.println(observationLines.size());
            }
        } catch (Exception e) {
            System.err.println(e);
        }
    }
}
