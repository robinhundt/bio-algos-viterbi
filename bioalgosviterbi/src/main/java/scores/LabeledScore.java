package scores;

public class LabeledScore {
    private final double score;
    private final boolean label;

    public LabeledScore(double score, boolean label) {
        this.score = score;
        this.label = label;
    }

    @Override
    public String toString() {
        return "Probability: " + String.valueOf(this.score) + "\tLabel: " + this.label;
    }

    /**
     * @return the score
     */
    public double getScore() {
        return score;
    }

    /**
     * @return the label
     */
    public boolean isLabel() {
        return label;
    }
}