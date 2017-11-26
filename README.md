ECPred Version 1.0 (01.12.2017)

Dependencies

Java 8 
gcc/g++ (any version)

Download

ECPred.jar

and

ECPred_training_data.tar.gz

Training data should be downloaded from:

..................

Installation

Run runLinux.sh or runMac.sh (while inside lib folder) from terminal using one of the following commands:

./runLinux.sh
or
./runMac.sh

These bash scripts will install necessary libraries and tools.

Usage

(ECPred.jar and lib folder should be at same directory)

Run java -jar ECPred.jar input.fasta

Input

ECPred accepts one input fasta file which may contain up to 20 proteins.

Output

ECPRed.jar outputs a tsv file predictionResults_inputFile_DateandTime.tsv which contains the main, subfamily, sub-subfamily and substrate class predictions together with confidence scores for each prediction; alternatively, the output can be “non-enzyme” or “no prediction” for each query protein, when there is no EC number prediction.

Data files

"Independent_test_results_1030_proteins.tsv": A tab delimited file containing uniProtIDs, true EC Number classes and ECPred predictions for the selected independent test set proteins (1030 sequences).

"Independent_test_results_60_proteins.tsv": A tab delimited file containing uniProtIDs, true EC Number classes and ECPred predictions for the selected independent test set proteins (60 sequences).

"ECNumberList.txt": A text file containing the list of EC numbers that ECPred can predict.

"test.fasta": An example input fasta file.

"predictionResults_test_DateandTime.tsv": An example output prediction file (for test.fasta).
