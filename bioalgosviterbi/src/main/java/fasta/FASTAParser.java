package fasta;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Optional;

public class FASTAParser {
    private static final int stringBuilderInitialCapacity = 80;

    public static ArrayList<FASTASequence> parse(Path path) throws IOException {
        var lines = Files.readAllLines(path);
        var sequences = new ArrayList<FASTASequence>();
        Optional<String> id = Optional.empty();
        StringBuilder stringBuilder = new StringBuilder(stringBuilderInitialCapacity);
        for (var i = 0; i < lines.size(); i++) {
            var currLine = lines.get(i);
            Optional<String> nextline = i+1 >= lines.size() ? Optional.empty() : Optional.of(lines.get(i+1));
            if (currLine.startsWith(">")) {
                id = Optional.of(currLine.substring(1)); // skip the >
                continue;
            }
            stringBuilder.append(currLine);
            if (!nextline.isPresent() || nextline.get().startsWith(">")) {
                if (!id.isPresent()) {
                    throw new IOException("Invalid FASTA File");
                }
                sequences.add(new FASTASequence(
                        id.get(), stringBuilder.toString().toCharArray()));
                stringBuilder = new StringBuilder(stringBuilderInitialCapacity);
            }
        }
        return sequences;
    }
}
