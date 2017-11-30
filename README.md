## ECPred Version 1.0 (01.12.2017)

## Dependencies

#### Java 8  

For Linux, you can install latest version of Java by running following commands from terminal:
```
sudo apt-get update
sudo apt-get install default-jre
sudo apt-get install default-jdk
```
For Mac, you can install latest version of Java by running following commands from terminal:
```
brew update
brew cask install java
```

#### g++ (any version)  

For Linux, you can install latest version of g++ by running following commands from terminal
```
sudo apt-get update
sudo apt-get install build-essential
```
For Mac, you can install latest version of g++ by running following commands from terminal  <br />

Type ``` g++ ``` on terminal.

## Download
```
ECPred.tar.gz
```
Above file (around 3 GB) should be downloaded from:

http://cansyl.metu.edu.tr/ECPred.html

## Installation

Extract the files using: <br />
```
tar -xvf ECPred.tar.gz  
```
After extraction the total size of the folder will be around 10 GB. <br />

Run runLinux.sh or runMac.sh from terminal according to your OS using one of these commands: <br />
```
./runLinux.sh 
```
or <br />
```
./runMac.sh
```
These bash scripts will install necessary libraries and tools.

## Usage

Run the following command on terminal to analyze the file "filename.fasta"  <br />
```
java -jar ECPred.jar filename.fasta
```
## Input

ECPred accepts one input fasta file which may contain up to 20 proteins.

## Output

"predictionResults_filename_Date-Time.tsv": <br />

A tsv file that contains the main, subfamily, sub-subfamily and substrate class predictions together with confidence scores for each prediction; alternatively, the output can be “non-enzyme” or “no prediction”.

## Data files

"ECNumberList.txt":  <br />

A text file containing the list of EC numbers that ECPred can predict.  <br />

"test.fasta":  <br />

An example input fasta file.  <br />

"predictionResults_input_20171128-172728.tsv":  <br />

An example output prediction file (for test.fasta).
