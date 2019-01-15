package app;

import fasta.FASTAParser;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
        try {

            var test = FASTAParser.parse(Paths.get("../data/LSU_short_f.fasta"));
            System.out.println(test);
            System.out.println(test.size());
        } catch (Exception e) {
            System.err.println(e);
        }
    }
}
