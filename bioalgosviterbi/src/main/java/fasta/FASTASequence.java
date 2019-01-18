package fasta;

import java.util.Arrays;
import java.util.Map;

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

    public int[] parseBasesToInt(Map<Character, Integer> map) {
        int length = this.sequence.length;
        int[] intBases = new int[length];
        for (var i = 0; i < length; i++) {
            intBases[i] = map.get(this.sequence[i]);
        }
        return intBases;
    }
}
