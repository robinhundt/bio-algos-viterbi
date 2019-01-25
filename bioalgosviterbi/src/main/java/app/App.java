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
            HashMap<Character, Integer> observationMap = new HashMap<Character, Integer>(Map.of('A', 0, 'C', 1, 'G', 2, 'U', 3));
            // import training sequences
            var trainingSequences = FASTAParser.parse(Paths.get(parameter.getTraining()));

            // create profil HMM with test sequences
            ProfileHMM pHMM = new ProfileHMM(trainingSequences, gapSymbol, observationMap, parameter.getEmissionPseudocounts()/*,parameter.getTransitionPseudocounts()*/,0.5);
            // double[][] emissionMatrix = pHMM.getEmissionMatrix();
            // double[][] transitionMatrix = pHMM.getTransitionMatrix();

            // System.out.println(printMatrix(emissionMatrix));
            // System.out.println(printMatrix(transitionMatrix));
            //get test sequences
            ArrayList<String> testFiles = getFileList(parameter.getTest());

            runVitberiOnTestFiles(parameter, observationMap, testFiles, pHMM);

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
            
            Scanner reader = new Scanner(System.in);
            PriorityQueue<LabeledScore> allLabeledScores = new PriorityQueue<LabeledScore>();
            
            int nrPositives = 0;
            int nrNegatives = 0;
            
            int countTruePositives = 0;
            int countFalsePositives = 0;
            ArrayList<Double> truePositiveRate = new ArrayList<Double>();
            ArrayList<Double> falsePositiveRate = new ArrayList<Double>();
            double fPrev = Double.MAX_VALUE;
            
            // get scores and labels for each file and store them in ascending order
            for (var filename : resultFiles) {
                System.out.println("Label for file:\n" + filename);
                var label = Boolean.parseBoolean(reader.next());

                PriorityQueue<LabeledScore> labeledScore = getLabeledScores(filename, label);
                
                if (label) {
                    nrPositives += labeledScore.size();
                } else {
                    nrNegatives += labeledScore.size();
                }
                allLabeledScores.addAll(labeledScore);
            }
            
            
            while (!allLabeledScores.isEmpty()) {
                var score = allLabeledScores.poll();
                if (score.getScore() != fPrev) {
                    truePositiveRate.add(countTruePositives/((double) nrPositives));
                    falsePositiveRate.add(countFalsePositives/((double) nrNegatives));
                    fPrev = score.getScore();
                }

                if (score.isLabel()) {
                    countTruePositives++;
                }
                else {
                    countFalsePositives++;
                }
            }
            System.out.println(truePositiveRate);

            // write TPR and FPR to file and create roc curve and AUC value with python
            List<String> lines = Arrays.asList(truePositiveRate.toString(), falsePositiveRate.toString());
            File rates = new File(parameter.getOutputFolder() + "rocCurve/" + "rocCurveData.txt");
            rates.createNewFile();
            Files.write(Paths.get(rates.getAbsolutePath()), lines, Charset.forName("UTF-8"));

            reader.close();

        } catch (Exception e) {
            System.out.println(e);
        }
    }

    private static PriorityQueue<LabeledScore> getLabeledScores(String resultFile, boolean label)
            throws IOException {

        PriorityQueue<LabeledScore> labeledScores = new PriorityQueue<LabeledScore>();
        Path path = Paths.get(resultFile);
        
        List<String> lines = Files.readAllLines(path);
        for (var line : lines) {
            // line looks like: viterbiPath;maxProbability
            LabeledScore score = new LabeledScore(Double.parseDouble(line.split(";")[1]), label);
            labeledScores.add(score);
        }

        return labeledScores;
    }

    private static void runVitberiOnTestFiles(Parameter parameter, HashMap<Character, Integer> observationMap,
            ArrayList<String> testFiles, ProfileHMM pHMM) throws IOException {

        for (var testFile : testFiles) {
            var sequences = FASTAParser.parse(Paths.get(testFile));
            ArrayList<ViterbiResult> vitProbabilities = new ArrayList<ViterbiResult>(testFile.length());

            // Random random = new Random();
            for (var sequence : sequences) {
                int[] observations = sequence.parseBasesToInt(observationMap);
          
                // calculate viterbi path and probability
                
                ViterbiResult viterbiResult = Viterbi.calc(observations, pHMM);
                // int[] path = {1,2,0,1,0,0,1,1,0};
                // double maxProbability = random.nextDouble();
                vitProbabilities.add(viterbiResult);
            }
            // create new file and store each probability in one line

            List<String> lines = new ArrayList<String>();
            for (var vitProbability : vitProbabilities) {
                lines.add(vitProbability.toString());
            }
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
