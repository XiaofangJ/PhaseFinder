# PhaseFinder

## Overview
The PhaseFinder algorithm is designed to detect DNA inversion mediated phase variation in bacterial genomes using genomic or metagenomic sequencing data. It works by identifying regions flanked by inverted repeats, mimicking their inversion in silico, and identifying regions where sequencing reads support both orientations. Here, we define phase variation as "a process employed by bacteria to generate frequent and reversible changes within specific hypermutable loci, introducing phenotypic diversity into clonal populations”. Not every region detected by PhaseFinder will directly result in phase variation, but the results should be highly enriched for regions that do. 

## Prerequisites
+ [Biopython](https://biopython.org/)
+ [pandas](https://pandas.pydata.org)
+ [samtools](http://samtools.sourceforge.net/) (>=1.4)
+ [bowtie](https://github.com/BenLangmead/bowtie)(>=version 1.2.0)
+ [einverted](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html)
+ [bedops](https://bedops.readthedocs.io/en/latest/)
+ [bedtools](https://bedtools.readthedocs.io/en/latest/)

On Ubuntu 18.04 the required software can be installed with the following commands.
```
sudo apt-get install bowtie samtools bedtools emboss python-pip python-pandas
sudo pip install biopython 
wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 && tar jxf bedops_linux_x86_64-v2.4.35.tar.bz2 && sudo cp bin/* /usr/local/bin/
```

## Quick Start
All you need to get started is a genome (in fasta format) you would like to search for invertible DNA regions and genomic sequencing data (preferrably Illumina in fastq format) from the same organism, or metagenomic sequencing data from a sample containing the organism (preferrably Illumina in fastq format). 

To test PhaseFinder, you can use the example files (genome: test.fa, genomic data: p1.fq, p2.fq) 

Example:
```
# Identify regions flanked by inverted repeats 
python PhaseFinder.py locate -f ./data/test.fa -t ./data/test.einverted.tab -g 15 85 -p 

# Mimic inversion
python PhaseFinder.py create -f ./data/test.fa -t ./data/test.einverted.tab -s 1000 -i ./data/test.ID.fasta

# Identify regions where sequencing reads support both orientations 
python PhaseFinder.py ratio -i ./data/test.ID.fasta -1 ./data/p1.fq -2 ./data/p2.fq -p 16 -o ./data/out
```

If successful, the output will be in data/out.ratio.txt

In this example, there is one real example of an invertible DNA region "am_0171_0068_d5_0006:81079-81105-81368-81394" because only this region has reads supporting both the F and R orientation. 

---

## Tutorial
### 1. Generate a position table of regions flanked by inverted repeats 
Users can identify inverted repeats using the "PhaseFinder.py locate" command, or generate their own table.

#### 1.1. Generate the position table with the PhaseFinder script
```
usage: PhaseFinder.py locate [-h] -f  -t  [-e] [-m] [-r] [-g ] [-p]

optional arguments:
  -h, --help         show this help message and exit
  -f , --fasta       input genome sequence file in fasta format
  -t , --tab         output table with inverted repeats coordinates
  -e , --einv        einverted parameters, if unspecified run with PhaseFinder
                     default pipeline
  -m , --mismatch    max number of mismatches allowed between IR pairs,used
                     with -einv (default:3)
  -r , --IRsize      max size of the inverted repeats, used with -einv
                     (default:50)
  -g  , --gcRatio    the minimum and maximum value of GC ratio
  -p, --polymer      Remove homopolymer inverted repeats
```

##### Input: A fasta file containing the genome sequence
##### Output: A table file containing the postion information of invereted repeats in the genome

##### Examples:
* Run the default PhaseFinder locate parameters
```
python PhaseFinder.py locate -f ./data/test.fa -t ./data/test.einverted.tab 
```
* Run the default PhaseFinder locate parameters and remove inverted repeats with GC content lower than 15% and higher than 85% or with homopolymers
```
python PhaseFinder.py locate -f ./data/test.fa -t ./data/test.einverted.tab -g 15 85 -p 
```
* Run with the specified einverted parameters "-maxrepeat 750 -gap 100 -threshold 51 -match 5 -mismatch -9" 
```
python PhaseFinder.py locate -f ./data/test.fa -t ./data/test.einverted.tab -e "-maxrepeat 750 -gap 100 -threshold 51 -match 5 -mismatch -9" 
```


#### 1.2. Generate the position table with other tools
You can identify regions flanked by inverted repeats directly with tools such as [einverted](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/einverted.html) and [palindrome](http://emboss.sourceforge.net/apps/cvs/emboss/apps/palindrome.html). 

![inverted repeats](https://github.com/XiaofangJ/PhaseFinder/blob/master/IR.png)

Prepare the output into the following format:

A table file with five columns (tab delimited):

 Column name | Explanation                                                   |
-------------|---------------------------------------------------------------|
 Scaffold    | The scaffold or contig name where the inverted repeat is detected
 pos A       | The start coordinate of the first inverted repeat (0-based)
 pos B       | The end coordinate of the first inverted repeat (1-based)
 pos C       | The start coordinate of the second inverted repeat (0-based)
 pos D       | The end coordinate of the second inverted repeat (1-based)

---
### 2. Mimic inversion in silico to create a database of inverted sequences
```
usage: PhaseFinder.py create [-h] -f  -t  -s  -i

optional arguments:
  -h, --help         show this help message and exit
  -f , --fasta       input genome sequence file in fasta format
  -t , --tab         table with inverted repeat coordinates
  -s , --flanksize   base pairs of flanking DNA on both sides of the
                     identified inverted repeats
  -i , --inv         output path of the inverted fasta file
```

#### Input
* The position table from step 1

#### Output
* A fasta file containing inverted (R) and non-inverted (F) putative invertible DNA regions flanked by sequences of specified length (bowtie indexed)
* A table file (with suffix ".info.tab") describing the location of inverted repeats in the above fasta file

---
### 3. Align sequence reads to inverted sequence database and calculate the ratio of reads aligning to the F or R orienation. 
```
usage: PhaseFinder.py ratio [-h] -i  -1  -2  [-p] -o

optional arguments:
  -h, --help       show this help message and exit
  -i , --inv       input path of the inverted fasta file
  -1 , --fastq1    first pair in fastq
  -2 , --fastq2    second pair in fastq
  -p , --threads   number of threads
  -o , --output    output prefix
```

#### Input
* Output from step 2
* fastq file of genomic or metagenomic sequence used to verify DNA inversion
* Number of threads used for bowtie alignment and samtools process
#### Output
* A table file (with suffix ".ratio.txt") containing the reads that supporting either R or F orientation of invertible DNA

 Column name | Explanation                                                                 |
-------------|-----------------------------------------------------------------------------|
Sequence     | Putative invertible regions(Format:Scaffold:posA-posB-posC-posD)
Pe_F         | The number of reads supprting the F orientation with paired-end information
Pe_R         | The number of reads supprting the R orientation with paired-end information
Pe_ratio     | Pe_R/(Pe_F + Pe_R). The percent of reads supporting the R orientation with the paired-end method
Span_F       | The number of reads supporting the F orientation spanning the inverted repeat by at least 10 bp on either side
Span_R       | The number of reads supporting the R orientation spanning the inverted repeat by at least 10 bp on either side
Span_ratio   | Span_R/(Span_F + Span_R). The percent of reads supporting the R orientation with the spanning method. 

True invertible regions have reads supporting both the F and R orientation. We recommend combining the information from both the paired-end (Pe) and spanning (Span) methods to find valid invertible DNA regions. Our default is to classify a region as invertible if at least 1% of reads support the R orientation with a minimum Pe_R > 5 and Span_R > 3. 

### 4. (Optional) Subset for intergenic invertible DNA regions 

If you are especially interested in intergenic regulatory regions, such as promoters, you can remove predicted invertible regions overlapping with coding sequences (CDS). First, obtain an annotation for the genome of interest from the NCBI or that you genereate yourself in GFF3 format. Second, subsubset the annotation for CDS regions only. Third, use the following command to process the output of PhaseFinder step 3 to obtain a list of intergenic putative invertible DNA regions.

```
sed '1d' output_from_phasefinder.ratio.txt| awk '{print $1"\t"$0}'|sed 's/:/\t/;s/-[^\t]*-/\t/'|sortBed |closestBed  -a - -b annotation.gff  -d |awk '$20!=0{print $3}' > intergenic_IDR.txt
```

## Citation

