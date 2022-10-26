# LepMapWrapper
A wrapper for LepMap3 (Pasi Rastas https://doi.org/10.1093/bioinformatics/btx494).
Intended to ease the learning curve with lincage analysis using LepMap3
Basesd on (and borrowing from) the clairemerot/lepmap3_pipeline on github
## Introduction
Documentation for LepMap3 is not easy to find and there is a considerable learning curve.
This wrapper leads you through the process in an interactive manner and generates a whole
slew of output organized by timestamps so analyses can be repeated and compared.
It is still wery much a work in progress.  At the moment it only handles one family gracefully 
although it should be possible to force multi-family data trhough it.
There is still plenty of options to LepMap3 that are not accounted for so there is a lot of
room for improvement.

## Usage instructions
- Git clone LepMapWrapper 
- Get the dependencies sorted
- Prepare input files (see below)
The wrapper sets up three directories
* 00_scripts/
* 01_input_files/
* 02_genome_files/
It is convenient to copy or symlink your input files and indexed genome / stacks_catalog.fasta.gz into the 
respective directories (althoug it is not required).
Run like so:

> 00_scripts/00_LepMað-Wrapper.sh /path/to/filtered.snps.vcf /path/to/pedigree_file.tsv

The scripts should now lead you through the rest

## Dependencies
- Linux or MacOS
- Perl / Python 2.7 / and R > v4 (plus 
- samtools (some fairly recent release)
