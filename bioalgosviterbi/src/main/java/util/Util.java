package util;

import java.util.HashMap;

public class Util {
    /**
     * Transforms the provided matrix by an element wise natural logarithm in place.
     *
     * @param matrix 2d matrix that should be transformed.
     * @return the transformed matrix.
     */
    public static double[][] toLog(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = Math.log(matrix[i][j]);
            }
        }
        return matrix;
    }

    public static HashMap<Character, Integer> createObersavtionMap() {
        var observationMap = new HashMap<Character, Integer>();
        observationMap.put('A', 0);
        observationMap.put('C', 1);
        observationMap.put('G', 2);
        observationMap.put('T', 3);
        return observationMap;
    }
}
