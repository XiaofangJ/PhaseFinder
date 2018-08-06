# PhaseFinder

## Overview
PhaseFinder is a pipeline to identify and quantify DNA inversion in bacterial genome using genomic/metagenomic sequencing data.

## Software requirement
+ [samtools](http://samtools.sourceforge.net/) (>=1.4)
+ [bowtie](https://github.com/BenLangmead/bowtie)(>=version 4.8.0)

## Data preparation

![inverted repeats](https://github.com/XiaofangJ/PhaseFinder/blob/master/IR.png)
Identify regions flanked by inverted repeats with tools such as [einverted](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html) and [palindrome](http://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html).


Prepare the output into the following format:

a table file with five columns(separated by tab):

 Column name | Explanation                                                   |
-------------|---------------------------------------------------------------|
 Scaffolds   | The scaffolds where the inverted repeats are detected
 pos A       | The start coordinates of the first inverted repeat (0-based)
 pos B       | The end coordinates of the first inverted repeat (1-based)
 pos C       | The start coordinates of the second inverted repeat (0-based)
 pos D       | The end coordinates of the second inverted repeat (1-based)

## Quick Start
Example:
```
python PhaseFinder.py create -f ./data/test.fa -t ./data/test.einverted.tab -s 1000 -i ./data/test.ID.fasta
python PhaseFinder.py ratio -i ./data/test.ID.fasta -1 ./data/p1.fq -2 ./data/p2.fq -p 16 -o ./data/out
```

## Tutorial
### 1. Create a database of inverted sequences
```
usage: PhaseFinder.py create [-h] -f  -t  -s  -i

optional arguments:
  -h, --help         show this help message and exit
  -f , --fasta       input genome sequence file in fasta format
  -t , --tab         table with inverted repeats coordinates
  -s , --flanksize   flanked base pair
  -i , --inv         output path of the inverted fasta file
```
#### Input
* The genome sequence where region flanked by inverted repeats identified
* The table with inverted repeats information (as described in the data preparation)
* The base pairs sequence flanked the putative invertible regions (default:1000)

#### Output
* Fasta file containing inverted(R) and non-inverted(F) putative invertible DNA regions flanked by sequences of specified length

### 2. Align sequence reads to inverted sequence database and caculate the ratio
```
usage: PhaseFinder.py ratio [-h] -i  -1  -2  [-p] -o

optional arguments:
  -h, --help       show this help message and exit
  -i , --inv       input path of the inverted fasta file
  -1 , --fastq1    first pair in fastq
  -2 , --fastq2    second pair in fastq
  -p , --threads   number of threads used
  -o , --output    output prefix
```
#### Input
* Output from step 1
* fastq file used to verify DNA inversion
* Number of threads used for bowtie alignment and samtools process
#### Output
* A table file (with suffix ".ratio.txt") containing the reads that supporting either R or F orientation of invertible DNA

 Column name | Explanation                                                                 |
-------------|-----------------------------------------------------------------------------|
Sequence     | Putative invertible regions(Format:Scaffold:posA-posB-posC-posD)
Pe_F         | The number of reads supprting F orientation with pair-end information
Pe_R         | The number of reads supprting R orientation with pair-end information
Span_F       | The number of reads supporting F orientation via spannign invertible region
Span_R       | The number of reads supporting R orientation via spannign invertible region


This table can be further process to select valid invertible region and  estimate the ratio of bacteria with R/F orientation in a population.

## Citation
