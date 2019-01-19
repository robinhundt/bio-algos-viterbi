package app;

import fasta.FASTAParser;
import parameter.Parameter;
import phmm.ProfileHMM;
import viterbi.*;
import scores.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.nio.charset.*;
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
            // import training sequences
            var trainingSequences = FASTAParser.parse(Paths.get(parameter.getTraining()));

            // create profil HMM with test sequences
            ProfileHMM pHMM = new ProfileHMM(trainingSequences, gapSymbol, observationMap, parameter.getEmissionPseudocounts()/*,parameter.getTransitionPseudocounts()*/);
            double[][] emissionMatrix = pHMM.getEmissionMatrix();
            double[][] transitionMatrix = pHMM.getTransitionMatrix();

            //get test sequences
            ArrayList<String> testFiles = getFileList(parameter.getTest());

            runVitberiOnTestFiles(parameter, observationMap, testFiles);

            // create roc curve
            if (parameter.isRocCurve()) {
                rocCurve(parameter);
            }
        } catch (Exception e) {
            System.err.println(e);
        }
    }

    private static void rocCurve(Parameter parameter) {
        ArrayList<String> resultFiles = getFileList(parameter.getOutputFolder());
        try {
            ArrayDeque<LabeledScore> labeledScores = getLabeledScores(resultFiles);
            int nrPositives = 0; // ?
            int nrNegatives = 0; // ?
            int countTruePositives = 0;
            int countFalsePositives = 0;
            ArrayList<Double> truePositiveRate = new ArrayList<Double>();
            ArrayList<Double> falsePositiveRate = new ArrayList<Double>();
            double fPrev = Double.MAX_VALUE;
            
            for (var score : labeledScores) {
                if (score.getScore() != fPrev) {
                    truePositiveRate.add(countTruePositives/(double) nrPositives);
                    falsePositiveRate.add(countFalsePositives/(double) nrNegatives);
                    fPrev = score.getScore();
                }

                if (score.isLabel()) {
                    countTruePositives++;
                }
                else {
                    countFalsePositives++;
                }
            }

            // write TPR and FPR to file and create roc curve and AUC value with python
            List<String> lines = Arrays.asList(truePositiveRate.toString(), falsePositiveRate.toString());
            Path path = Paths.get(parameter.getOutputFolder() + "rocCurve/" + "rocCurveData.txt");
            File rates = new File(path.toString());
            rates.createNewFile();
            Files.write(path, lines, Charset.forName("UTF-8"));

        } catch (Exception e) {
            System.out.println(e);
        }
    }

    private static ArrayDeque<LabeledScore> getLabeledScores(ArrayList<String> resultFiles)
            throws IOException {
        ArrayDeque<LabeledScore> allLabeledScores = new ArrayDeque<LabeledScore>();
        for (var filename : resultFiles) {
            // ask user for label
            ArrayDeque<LabeledScore> labeledScores = new ArrayDeque<LabeledScore>();
            // HashMap<Double, Boolean> labeledScores = new HashMap<Double, Boolean>();
            Path path = Paths.get(filename);
            System.out.println(filename);
            List<String> lines = Files.readAllLines(path);
            for (var line : lines) {
                LabeledScore score = new LabeledScore(Double.parseDouble(line), true);
                labeledScores.add(score);
            }
            allLabeledScores.addAll(labeledScores);
        }
        return allLabeledScores;
    }

    private static void runVitberiOnTestFiles(Parameter parameter, HashMap<Character, Integer> observationMap,
            ArrayList<String> testFiles) throws IOException {

        for (var testFile : testFiles) {
            var sequences = FASTAParser.parse(Paths.get(testFile));
            ArrayList<int[]> observationLines = new ArrayList<int[]>(sequences.size());
            ArrayList<String> vitProbabilities = new ArrayList<String>(testFile.length());

            for (var sequence : sequences) {
                observationLines.add(sequence.parseBasesToInt(observationMap));
            }

            for (var observations : observationLines) {
                // calculate viterbi path and probability
                /*
                 * ViterbiResult viterbiResult = Viterbi.calc(observations, transitionMatrix,
                 * emissionMatrix);
                 * vitProbabilities.add(String.valueOf(viterbiResult.getMaxProbability()));
                 */
                ;
            }
            // create new file and store each probability in one line

            List<String> lines =  vitProbabilities;
            String[] file = testFile.split("/");
            String filename = file[file.length - 1];
            File results = new File(parameter.getOutputFolder() + "probabilities-" + filename);
            results.createNewFile();
            Files.write(Paths.get(results.getAbsolutePath()), lines, Charset.forName("UTF-8"));

            // System.out.println(observationLines.get(0));
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
        boolean rocCurve = false;
        for (var line: lines) {
            if (line.startsWith("//")) {
                continue;
            }
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
                case "rocCurve": rocCurve = Boolean.parseBoolean(setting[1]); break;
            }
        }
        return new Parameter(training, test, outputFolder, backgroundFile, emissionPseudocounts, transitionPseudocounts, rocCurve);
    }

    private static ArrayList<String> getFileList(String path) {
        ArrayList<String> testFiles = new ArrayList<String>();
        Path filePath = Paths.get(path);

        if (Files.isDirectory(filePath)) {
            File directory = new File(path);

            for (var file : directory.listFiles()) {
                if (file.isFile()) {
                    testFiles.add(file.getAbsolutePath());
                }
            }
        } else if (Files.isRegularFile(filePath)) {
            testFiles.add(path);
        }

        return testFiles;
    }
}
