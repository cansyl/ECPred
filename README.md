## ECPred Version 1.0 (01.12.2017)

## Dependencies

Java 8 <br />
gcc/g++ (any version)

## Download

ECPred.jar <br />

and <br />

ECPred_training_data.tar.gz <br />

Training data should be downloaded from: <br />

..................

## Installation

Run runLinux.sh or runMac.sh (while inside lib folder) from terminal using one of the following commands: <br />

./runLinux.sh <br />
or <br />
./runMac.sh <br />

These bash scripts will install necessary libraries and tools.

## Usage

(ECPred.jar and lib folder should be at same directory) <br />

Run java -jar ECPred.jar input.fasta

## Input

ECPred accepts one input fasta file which may contain up to 20 proteins.

## Output

ECPRed.jar outputs a tsv file predictionResults_inputFile_Date-Time.tsv which contains the main, subfamily, sub-subfamily and substrate class predictions together with confidence scores for each prediction; alternatively, the output can be “non-enzyme” or “no prediction” for each query protein, when there is no EC number prediction.

## Data files

"Independent_test_results_1030_proteins.tsv": A tab delimited file containing uniProtIDs, true EC Number classes and ECPred predictions for the selected independent test set proteins (1030 sequences). <br />

"Independent_test_results_60_proteins.tsv": A tab delimited file containing uniProtIDs, true EC Number classes and ECPred predictions for the selected independent test set proteins (60 sequences). <br />

"ECNumberList.txt": A text file containing the list of EC numbers that ECPred can predict. <br />

"test.fasta": An example input fasta file. <br />

"predictionResults_test_20171127-012847.tsv": An example output prediction file (for test.fasta).
