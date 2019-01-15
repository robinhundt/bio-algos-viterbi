package fasta;

public class FASTASequence {
    private final String id;
    private final char[] sequence;

    public FASTASequence(String id, char[] sequence) {
        this.id = id;
        this.sequence = sequence;
    }

    @Override
    public String toString() {
        return String.valueOf(sequence);
    }

    public String getId() {
        return id;
    }

    public char[] getSequence() {
        return sequence;
    }
}
