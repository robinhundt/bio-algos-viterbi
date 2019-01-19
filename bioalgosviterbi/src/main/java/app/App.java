package app;

import fasta.FASTAParser;
import parameter.Parameter;
import phmm.ProfileHMM;
import viterbi.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        if (args.length < 1) {
            System.out.println("Helpmessage");
            return;
        }
        try {
            Parameter parameter = parseParameters(args[0]);

            Character gapSymbol = '-';
            HashMap<Character, Integer> observationMap = new HashMap<Character, Integer>(Map.of('A', 0, 'C', 1, 'G', 2, 'U', 3, gapSymbol, 4));
            int pseudoCount = 1;
            // import training sequences
            var trainingSequences = FASTAParser.parse(Paths.get(parameter.getTraining()));

            // create profil HMM with test sequences
            ProfileHMM pHMM = new ProfileHMM(trainingSequences, gapSymbol, observationMap, pseudoCount);
            double[][] emissionMatrix = pHMM.getEmissionMatrix();
            double[][] transitionMatrix = pHMM.getTransitionMatrix();
            for (var column : emissionMatrix) {
                System.out.println(Arrays.toString(column));
                
            }

            //get test sequences
            ArrayList<String> testFiles = getFileList(parameter.getTest());

            for (var testFile : testFiles) {
                var sequences = FASTAParser.parse(Paths.get(testFile));
                ArrayList<int[]> observationLines = new ArrayList<int[]>(sequences.size());
                ArrayList<String> vitProbabilities = new ArrayList<String>(testFile.length());

                for (var sequence : sequences) {
                    observationLines.add(sequence.parseBasesToInt(observationMap));
                }

                for (var observations: observationLines) {
                    //calculate viterbi path and probability
                    /*
                    ViterbiResult viterbiResult = Viterbi.calc(observations, transitionMatrix, emissionMatrix);
                    vitProbabilities.add(String.valueOf(viterbiResult.getMaxProbability()));
                    */
                    ;
                }
                // create new file and store each probability in one line
                String[] file = testFile.split("/");
                String filename = file[file.length-1];
                File results = new File(parameter.getOutputFolder() + "probabilities-" + filename);
                results.createNewFile();
                //Files.write(results.getAbsolutePath(), vitProbabilities);;
                

                // System.out.println(observationLines.get(0));
            }
        } catch (Exception e) {
            System.err.println(e);
        }
    }

    private static Parameter parseParameters(String parameterFile) throws IOException {
        List<String> lines = Files.readAllLines(Paths.get(parameterFile));
        String training = "";
        String test = "";
        String outputFolder = "";
        String backgroundFile = "";
        int emissionPseudocounts = 1;
        int transitionPseudocounts = 1;
        for (var line: lines) {
            var setting = line.split(":");
            if (setting.length < 2) {
                continue;
            }
            switch (setting[0]) {
                case "trainingData": training = setting[1]; break;
                case "testData": test = setting[1]; break;
                case "outputFolder": outputFolder = setting[1].endsWith("/") ? setting[1] : setting[1] + "/"; break;
                case "backgroundFile": backgroundFile = setting[1]; break;
                case "emissionPseudocounts": emissionPseudocounts = Integer.parseInt(setting[1]); break;
                case "transitionPseudocounts": transitionPseudocounts = Integer.parseInt(setting[1]); break;
            }
        }
        return new Parameter(training, test, outputFolder, backgroundFile, emissionPseudocounts, transitionPseudocounts);
    }

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
}
