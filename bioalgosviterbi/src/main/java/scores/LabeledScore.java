package scores;

public class LabeledScore implements Comparable {
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

    @Override
    public int compareTo(Object o) {
        if (o instanceof LabeledScore) {
            int comp = 0;
            if (this.score > ((LabeledScore) o).score) {
                comp = 1;
            } else if (this.score < ((LabeledScore) o).score) {
                comp = -1;
            }
            return comp;
        } else {
            throw new ClassCastException();
        }
    }
}