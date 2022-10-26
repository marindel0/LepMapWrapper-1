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

>\> 00_scripts/00_LepMap-Wrapper.sh /path/to/filtered.snps.vcf /path/to/pedigree_file.tsv

The "master script" should now lead you through the rest.

## Dependencies
- Linux or MacOS
- Perl / Python 2.7 / and R > v4 (plus packages 
- samtools (some fairly recent release)
Different parts of the pipeline will depend on different things (this is a mix & match) so most of it should
run even if some dependencies are missing.  If things go wrong, just try again.

## Preparation of input files
We need two types of input files to get the pipeline started:

1) A vcf file with filtered snps with high genotyping coverage in the family, and 
2) A custom pedigree file that is a bit of manual labor to assemble (see below).

LepMap3 allows other types of input but vcf is what this wrapper is customized for.  Later stages involve
retrieving DNA sequence from either a genome or a stacks catalog.fasta.gz file from Stacks. To run them 
you will need either (or both):

3) The genome fasta file that was used for sequence alignment / snp calling, which can be bgzip compressed.
* The genome file should be faidx indexed if samtools are used for seq extraction
4) If snps were called using the Stacks denovo pipeline then sequences can be retrieved from the catalog.fasta.gz
* The catalog should already be gzipped - don't unzip it.



