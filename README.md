# ECPred Version 1.0 26/11/2017

## Dependencies
Java 8 <br />
Any version of gcc/g++

## Installation
Inside of lib folder, Run runLinux.sh or runMac.sh from terminal using one of the following commands, <br />

./runLinux.sh or ./runMac.sh. <br />

These bash scripts will install necessary libraries and tools.

## Usage
ECPred.jar and lib folder should be at same directory. <br /> 
Run java -jar ECPred.jar input.fasta

## Input
ECPred accepts one input fasta file which may contain up to 20 proteins.

## Output
ECPRed.jar outputs a tsv file predictionResults_inputFile_DateandTime.tsv which contains the main, <br />
subfamily, sub-subfamily and substrate class predictions together with confidence scores for each <br />
prediction; alternatively, the output can be “non-enzyme” or  “no prediction” for each query protein, <br />
when there is no EC number prediction.

## Independent Performance Test Results
"SwissProt_v2017_06_predictions.tsv": A tab delimited file containing uniProtIDs, true EC Number classes <br />
and ECPred predictions for 1030 proteins.

## Independent Performance Test Results2
"SwissProt_v2017_06_predictions2.tsv": A tab delimited file containing uniProtIDs, true EC Number classes <br />
and ECPred predictions for selected 60 proteins.

## EC Number List
"ECNumberList.txt": A text file containing EC numbers that ECPred predicts.
