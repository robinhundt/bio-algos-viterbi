# Implementation

Dieser Ordner enthält die Implementation der HMM und PHMM Strukturen sowie jeweils einen Viterbitalgorithmus.

## Requirements
mind. Java 10

## Build
Zum Bauen sowie Testen des Projektes wird Maven verwendet:
```
mvn package
```
führt die in  `src/test` liegenden Tests aus, kompiliert das Programm und packt es in ein `.jar` Archiv.  
Bei Ausführung des Archivs muss der Pfad zu einer Konfigurationsdatei übergeben werden, Bsp. in `../data/parameters.txt`.  
Der Aufruf sieht also so aus:
```
java -jar target/bioalgosviterbi-0.1.0.jar ../data/parameters.txt
```
Über die Werte in der Parameterdatei können folgende Hyperparameter angepasst werden:
- trainingData
- testData
- outputFolder
- emissionPseudocounts
- transitionPseudocounts
- deleteDeletePseudocounts
- rocCurve


## Anmerkung!
Auf Grund der hohen Laufzeitklasse des Viterbialgorithmus, `O(|Observations|*|States|)`, bzw. `O(|States|^2)` benötigt der Viterbialgorithmus für das ProfilHMM eine enorme Menge an Arbeitsspeicher (teilweise über 4GB). 