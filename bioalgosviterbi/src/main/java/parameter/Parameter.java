package parameter;

public class Parameter {
    private final String training;
    private final String test;
    private final String outputFolder;
    private final String backgroundFile;
    private final int emissionPseudocounts;
    private final int transitionPseudocounts;
    private final int deleteDeletePseudocounts;
    private final boolean rocCurve;

    public Parameter(String training, String test, String outputFolder, String backgroundFile, int emissionPseudocounts, int transitionPseudocounts, int deleteDeletePseudocounts, boolean rocCurve) {
        this.training = training;
        this.test = test;
        this.outputFolder = outputFolder;
        this.backgroundFile = backgroundFile;
        this.emissionPseudocounts = emissionPseudocounts;
        this.transitionPseudocounts = transitionPseudocounts;
        this.deleteDeletePseudocounts = deleteDeletePseudocounts;
        this.rocCurve = rocCurve;
    }

    /**
     * @return the training
     */
    public String getTraining() {
        return training;
    }

    /**
     * @return the test
     */
    public String getTest() {
        return test;
    }

    /**
     * @return the outputFolder
     */
    public String getOutputFolder() {
        return outputFolder;
    }

    /**
     * @return the backgroundFile
     */
    public String getBackgroundFile() {
        return backgroundFile;
    }

    /**
     * @return the emissionPseudocounts
     */
    public int getEmissionPseudocounts() {
        return emissionPseudocounts;
    }

    /**
     * @return the transitionPseudocounts
     */
    public int getTransitionPseudocounts() {
        return transitionPseudocounts;
    }

    /**
     * @return the rocCurve
     */
    public boolean isRocCurve() {
        return rocCurve;
    }

    /**
     * @return the deleteDeletePseudocounts
     */
    public int getDeleteDeletePseudocounts() {
        return deleteDeletePseudocounts;
    }

}