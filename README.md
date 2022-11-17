# LepMapWrapper
A wrapper for LepMap3 (Pasi Rastas https://doi.org/10.1093/bioinformatics/btx494).
Intended to ease the learning curve with linkage analysis using LepMap3
Basesd on (and borrowing heavily from) the [clairemerot/lepmap3_pipeline](https://github.com/clairemerot/lepmap3_pipeline) on github

## Introduction
There is a lot of valuable information in the [LepMap3 wiki](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/)
but there is a considerable learning curve in order to get things to work.
This wrapper tries to lead you through the process in an interactive manner and generates a whole
slew of output organized by timestamps so analyses can be repeated and compared, in a convenient way.
It is still wery much a work in progress.  At the moment it only handles one family gracefully
although it should be possible to force multi-family data trhough it.
There is still plenty of options to LepMap3 and modules that are not accounted for so there is a lot of
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
```console
$> 00_scripts/00_LepMap-Wrapper.sh /path/to/filtered.snps.vcf /path/to/pedigree_file.tsv
```
The "master script" should now lead you through the rest. You might have to make the script executable
`chmod +x 00_scripts/00_LepMap-Wrapper.sh` or simply add sh or bash in front.

## Dependencies
- [LepMap3](https://sourceforge.net/projects/lep-map3/) (obviously) and a working java runtime environment
- Linux or MacOS
- Perl 5+
- Python 2.7 (not strictly required)
- R > v4 (with the dplyr package installed)
- Samtools (a recent release)

Besides, vcftools are convenient for pre-filtering if following this guide.

Different parts of the pipeline will depend on different things (this is a mix & match) so most of it should
run even if some dependencies are missing.  If things go wrong, just try again.

## Preparation of input files
We need two types of input files to get the wrapper-pipeline started:

1) A vcf file containing filtered snps with high genotyping coverage in the family, and
2) A custom pedigree file that is a bit of manual labor to assemble (see below).

LepMap3 allows several types of input, but vcf is what this wrapper is customized for.  Later stages involve
retrieving DNA sequence from either a *genome* or a stacks `catalog.fasta.gz` file from Stacks. To run them
you will need either (or both):

3) The genome fasta file that was used for sequence alignment / snp calling, which can be bgzip compressed.
* The genome file should be faidx indexed if samtools are used for seq extraction
4) If snps were called using the Stacks denovo pipeline then sequences can be retrieved from the catalog.fasta.gz
* The catalog should already be gzipped - don't unzip it.

### The pedigree file
LepMap3 requires a custom pedigree file that can include one or more families with or without grandparents etc.
The format is "transposed" and **tab separated**. It is not easy to show the tabbed format in Readme file but to summarize:
- There will be six lines in the table.  All lines start with CHR \t POS \t (for some strange reason).
- There will be one column per individual (plus the two columns with CHR and POS)

* Line 1: Family: For example "1# for a single family but if multiple families are analyzed at once this will vary
* Line 2: Id (Individual): The name of the individual "sample"
* Line 3: Father:  The name of the father (but enter 0 in the columns for father and mother unless there are grandparents)
* Line 4: Mother:  The name of the mother (same as above)
* Line 5: Sex: Male = 1, Female = 2, Unknown = 0.  Can be 0 for all except the parents.
* Line 6: Phenotype: Irrelevant for us. Enter 0 for all.

It is important to get the parents and grandparents right, the rest is simple.  
Here is an examlpe:
```
CHR POS 1    1    1    1     1     1     1     1     1      etc...
CHR POS Amma Afi  Baba Diado Mama  Papa  F2-1  F2-2  F2-3
CHR POS 0    0    0    0     Diado Afi   Papa  Papa  Papa
CHR POS 0    0    0    0     Baba  Amma  Mama  Mama  Mama
CHR POS 2    1    2    1     2     1     0     0     0
CHR POS 0    0    0    0     0     0     0     0     0
```
It is probably easiest to use a spreadsheat to edit and \"transpose paste\" from a Stacks catalog file, Radiator Strata file
or extract the header from the vcf file used as input.  There is a little perl script found floating around on the
internet called `transposeTabDelimited.pl` in `00_scripts/utilities/` that can be helpful.
One way to start is to get the list of names from the Stacks catalog like this:
```console
> zcat catalog.calls | head -n 15 | grep '^#CHROM' | cut --complement -f2-9 > file.txt
```
The vcf file can also be used the same way skipping the zcat step.

### The vcf file
Stacks populations can output VCFv4.2 format files that are convenient as input. It is probably best to do relatively loose
filtering in populations and then use vcftools to filter rather aggressively.  It is most convenient to * *write single snp* *
in Stacks. We only want high coverage loci.
Here is an example of a pre-filtering step:
```console
> vcftools --vcf AdamsFam.vcf --min-alleles 2 --max-alleles 4 --max-missing 0.9 --mac 15 --remove-indels -c --recode > AdamsFam_filtered.vcf
```
Note that the --max-missing parameter is a bit counter-intuitive, 0.9 actually means that max allowed missingness is 0.1 or 10%


### Other tricks
If you are using snps from reads mapped onto a related genome as input you may want to remove markers placed outside scaffolds.
The details of how to do this depend on the source of the genome because of different naming conventions used (ensembl / ncbi)
but here is an example of how this can be done for ncbi genomes wit scaffold names starting in "NC_"

* Vcf-sort, bgzip and tabix index the vcf file from Stacks populations:
```console
> vcf-sort AdamsFam_filtered.vcf > AdamsFam_filt_sort.vcf
> bgzip AdamsFam_filt_sort.vcf
> tabix AdamsFam_filt_sort.vcf.vcf.gz
```
* Generate the list of scaffolds (chromosomes) and unplaced contigs:
```console
> zcat AdamsFam_filt_sort.vcf.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > AdamsFam_all_contigs.txt
> grep "NC_" AdamsFam_all_contigs.txt > AdamsFam_scaffolds_only.txt
```
* We need the list of chromosomes as a single space separated line with --chr before each entry.
```console
> echo "" $(cat AdamsFam_all_contigs.txt) |sed -e 's/ / --chr /g' -- > chr-contiglist.txt
```
This list can then be used to filter the vcf like so.
```console
> vcftools --vcf AdamsFam_filtered.vcf --recode $(cat chr-contiglist.txt) --out AdamsFam_filt_scaffolds.vcf
```
### Other ways to generate input files

Another way to prepare input data for LepMap3 and, if I understand the paper correctly, the preferred way is to parse reference aligned bam files
directly into posterior probabilities.  The clairmerot/lepmap3 pipeline instructions refer to using the awk scripts pileupParser2.awk and pileup2posterior.awk
from LepMap3 like this.
```console
> samtools mpileup -q 10 -Q 10 -s $(cat sorted_bams.txt)|awk -f pileupParser2.awk|awk -f pileup2posterior.awk|gzip >all_fam_post.gz
```
where `sorted_bams.txt` is a list of the sorted bam files to be used as input. The directory where the pipe is run must also contain a list of sample names
under the name `mapping.txt`.  Running those scripts for a large number of samples can take a while and according to the LepMap3 wiki have been replaced by a java
class `Pileup2Likelihoods`. The old scripts seem to have been removed although the wiki still refers to them but the have been repurposed in another git repo 
[icruz1989/IBDcalculation](https://github.com/icruz1989/IBDcalculation).

Here are some shortcuts to getting this done.  I am assuming that there is one bam file per sample and that they are named by sample name.  This is not absolutely required but I won't get into that here.  See the LepMap3 wiki.

* First sort the reference-aligned bam files
```console
> mkdir path_for_sorted_bams
> for file in path_to_bams/*.bam; do samtools sort -@ 30 $file > path_for_sorted_bams/$(basename ${file%\.bam}_sorted.bam); done
```
The -@ option specifies the number of threads to use, so set that according to your resources.
* Next prepare the two text files to use as parameters:

a) List of the sorted bam files to use as input
```console
> ls -1 sorted_LB_aligned_bam/ |tr "\n" "\t/" >sorted_bams.txt
```
b) List of samplenames (essentially the same list without the "sorted.bam" extension) - named mapping.txt
```console
> sed -e 's/_sorted.bam//g' sorted_bams.txt >mapping.txt
```
* Then use samtools to feed Pileup2Likelihood like this:
```console
> samtools mpileup -q 10 -Q 10 -s $(cat sorted_bams)|java -cp bin/ Pileup2Likelihoods|gzip >post.gz
```
The [LepMap3 wiki](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/) shows how this can be parallelized by running the pipe contig by contig as otherwise this can take hours.
