#!/usr/bin/bash

printf "\nThe end of the unified pipeline script for running LepMap3 \n" 

### Make sure genome path is passed as parameter
if [ "$#" -lt 2 ]; then
    printf "\nSomething is missing. You need to pass a genome fasta file (with fasta.fai index)) and the path to the LM-Run you want to analyze as parameters \n"
    printf "Usage: $0 path_to_genome.fasta LM-Run-directory/   \n\n" 1>&2
    exit 1
fi

GENOME=$1
if [ -f "$GENOME" ]; then
    printf "\n$GENOME file found.\n"
else 
    echo "file $GENOME does not exist. Please try again!"
    exit 1
fi

LMRUNDIR=$2   # TODO - make this more robust

INDIR="${LMRUNDIR}08_analyze_maps/01_maps/"
MAPCOMPDIR="${LMRUNDIR}08_analyze_maps/04_prepare_map_comp/"
CHROMODIR="${LMRUNDIR}08_analyze_maps/05_prepare_chromonomer/"

mkdir -p $MAPCOMPDIR
mkdir -p $CHROMODIR

#this is a short pipeline in itself,
#first use R to format the map into en entry file for the .py that extracts a sequence around a given position
#then run the .py to extract the sequence,
#then use R to format for mapcomp

ls $INDIR | while read i
	do 
	echo $i
	Rscript 00_scripts/Rscripts/format_for_mapcomp_1st_step.R $i $INDIR $MAPCOMPDIR 
	python 00_scripts/utilities/01_extract_snp_variants_with_flanking_claire.py $GENOME ${MAPCOMPDIR}$i".snplist" 100 ${MAPCOMPDIR}$i".seq"
	Rscript 00_scripts/Rscripts/format_for_mapcomp_3rd_step.R $i $INDIR $MAPCOMPDIR
	done
	
	#perl -pe 's/\\*//' "08_analyze_maps/04_prepare_map_comp/"$i".snplist" > "08_analyze_maps/04_prepare_map_comp/"$i".snplist" #remove stars on sex markers
	echo "4th step - combining two map files must be run manually"   #TODO - look into this
	
#this is a short pipeline in itself,
#first use R to format the map into en entry file for the .py that extracts a sequence around a given position
#then run the .py to extract the sequence,
#the 3rd one use R again on the original map file to format for chromocomer


ls $INDIR | while read i
	do 
	echo $i
	Rscript 00_scripts/Rscripts/format_for_chromonomer_1st_step.R $i $INDIR $CHROMODIR
	python 00_scripts/utilities/01_extract_snp_variants_with_flanking_claire_into_fasta.py $GENOME ${CHROMODIR}$i".snplist" 100 ${LMRUNDIR}$i".fasta"
	Rscript 00_scripts/Rscripts/format_for_chromonomer_3rd_step.R $i $INDIR $CHROMODIR
	done
	
	#TODO - Find logfile and append to it
	
	
	
	
