#!/usr/bin/bash

printf "\nThe end of the unified pipeline script for running LepMap3 \n" 

### Make sure genome path is passed as parameter
if [ "$#" -lt 2 ]; then
    printf "\nSomething is missing. We need both genome (in indexed fasta) and PEDIGREE input files \n"
    printf "Usage: $0 path_to_vcf_file path_to_pedigreefile \n\n" 1>&2
    exit 1
fi

#this is a short pipeline in itself,
#first use R to format the map into en entry file for the .py that extracts a sequence around a given position
#then run the .py to extract the sequence,
#then use R to format for mapcomp


ls 08_analyze_maps/01_maps/ | while read i
	do 
	echo $i
	Rscript 01_scripts/Rscripts/format_for_mapcomp_1st_step.R "$i"
	python 01_scripts/utilities/01_extract_snp_variants_with_flanking_claire.py ../mapcomp/02_data/genome/genome.fasta "08_analyze_maps/04_prepare_map_comp/"$i".snplist" 100 "08_analyze_maps/04_prepare_map_comp/"$i".seq"
	Rscript 01_scripts/Rscripts/format_for_mapcomp_3rd_step.R "$i"
	done
	
	#perl -pe 's/\\*//' "08_analyze_maps/04_prepare_map_comp/"$i".snplist" > "08_analyze_maps/04_prepare_map_comp/"$i".snplist" #remove stars on sex markers
	
	
#this is a short pipeline in itself,
#first use R to format the map into en entry file for the .py that extracts a sequence around a given position
#then run the .py to extract the sequence,
#the 3rd one use R again on the original map file to format for chromocomer


ls 08_analyze_maps/01_maps/ | while read i
	do 
	echo $i
	Rscript 01_scripts/Rscripts/format_for_chromonomer_1st_step.R "$i"
	python 01_scripts/utilities/01_extract_snp_variants_with_flanking_claire_into_fasta.py ../mapcomp/02_data/genome/genome.fasta "08_analyze_maps/05_prepare_chromonomer/"$i".snplist" 100 "08_analyze_maps/05_prepare_chromonomer/"$i".fasta"
	Rscript 01_scripts/Rscripts/format_for_chromonomer_3rd_step.R "$i"
	done
	
	
	
	
