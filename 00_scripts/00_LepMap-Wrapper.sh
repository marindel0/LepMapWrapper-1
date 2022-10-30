#!/usr/bin/bash

### Zophonías O. Jónsson - zjons@hi.is - 15.Oct.2022
### Unsing ideas and code from github: clairemerot/lepmap3_pipeline 
RUN_STARTS=$(date '+%y%m%d@%H%M')

printf "\nLepMapWrapper: The unified pipeline script for running LepMap3 - starting at ${RUN_STARTS} \n" 


### Make sure input files are passed as parameters
if [ "$#" -lt 2 ]; then
    printf "\nSomething is missing. We need both GENOTYPES (in vcf format) and PEDIGREE input files \n"
    printf "Usage: $0 path_to_vcf_file path_to_pedigreefile \n\n" 1>&2
    exit 1
fi
INPUT=$1
PEDIGREE=$2
#echo "Inputfile: $INPUT"
if [ -f "$INPUT" ]; then
    printf "\n$INPUT exists.\n"
else 
    echo "file $INPUT does not exist. Please try again!"
    exit 1
fi
if [ -f "$PEDIGREE" ]; then
    printf "$PEDIGREE exists.\n"
    NUM_SAMPLES=$(($(awk -F'\t' '{print NF; exit}' $PEDIGREE)-2))
    printf "$NUM_SAMPLES in pedigree fil\n"
else 
    echo "file $PEDIGREE does not exist. Please try again!"
    exit 1
fi

### prepare paths and filenames
echo
BASENAME=$(basename $INPUT)
echo "Basename: $BASENAME"
#DIR=$(dirname ${INPUT})
echo "Relative path: $(dirname ${INPUT})"
#CORENAME="$(basename $INPUT .vcf)"
CORENAME=${BASENAME%%.*}
echo "Resulting files will include \"$CORENAME\" in their names."
OUTDIR="LM-Run-${RUN_STARTS}"
echo "Results and log will be written to ${OUTDIR}"

#Fix stupid behaviour of zcat on MacOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    function zcat { 
        command zcat < "$1" 
    }
fi

#make outdir - go there and start logging.
mkdir $OUTDIR

LOGFILE="${OUTDIR}/${RUN_STARTS}-Run_LepMap3-${CORENAME}.log"
echo "logfile (parameters): ${LOGFILE}"

{  #Set up the logfile
    printf "Logfile for LepMap3 pipeline run ${RUN_STARTS}\n\n" 
    printf "%-25s %s\n" "Parameter 1 (VCF file): " $INPUT
    printf "%-25s %s\n" "Parameter 2 (Pedigree): " $PEDIGREE
    printf "%-25s %s\n" "Working directory: " $(pwd)
    printf "%-25s %s\n" "Full path to VCF file: " $(dirname $(readlink -f $INPUT))
    printf "%-25s %s\n" "Basename of VCF file: " $BASENAME
    printf "%-25s %s" "Corename of VCF file: " $CORENAME
    printf "    (this will be used to name output files from various steps)\n\n"
}>$LOGFILE

### Input file is there, now prompt for more parameters or set defaults

printf "\nNow you will be asked to enter various parameters\n\n"
#echo "What was the minimum coverage set in previous filtering steps? This will be used for logfile only"
#read -p "Enter min_cov - default [10]: " MIN_COV
#MIN_COV=${MIN_COV:-10}
#echo $MIN_COV
#echo
echo "This is just for the logfile ..."
read -p "Enter name of family - default [TheSimpsons]: " FAMILY
FAMILY=${FAMILY:-TheSimpsons}
echo $FAMILY
echo

### How many cpu's should we use?
function get_cpu_count {
    #CPU_COUNT=$(nproc --all)
    CPU_COUNT=$(getconf _NPROCESSORS_ONLN)                       #make it work cross-platform 
    echo "Number of CPU threads detected on system = $CPU_COUNT"
    One_less=$(($CPU_COUNT-1))
    read -p "How many threads should we use? Suggested - [$One_less]: " CPU
    CPU=${CPU:-$One_less}
    echo "CPU set $CPU"
    if [ ${CPU} -gt ${CPU_COUNT} ]; then
        echo "That makes no sense. We will not use more than the number of available threads"
        CPU=$CPU_COUNT
    fi
    printf "$CPU processors will be used\n\n"
}
get_cpu_count


### Make sure that LepMap3 is present and accessible
read -p "Where is LepMap3 located? - default [/programs/LepMap3/bin/]: " LEPMAPDIR
LEPMAPDIR=${LEPMAPDIR:-/programs/LepMap3/bin/}
if [ -f "$LEPMAPDIR/ParentCall2.class" ]; then
    printf "Lepmap found - all is good \n"
else
    until [ -f "$LEPMAPDIR/ParentCall2.class" ]
    do
        read -p "LepMap3 not found.  Try again: " LEPMAPDIR
    done
    printf "Lepmap found - finally \n"
fi

{
    printf "%-25s %s %s %s %s\n" "CPU core use:" $CPU " out of" $CPU_COUNT "available"
    printf "%-25s %s\n\n" "Path to LepMap3:" $LEPMAPDIR
   # printf "%-25s %s\n" "Input min coverage:" $MIN_COV
    printf "%-25s %s\n" "Family name:" $FAMILY
}>>$LOGFILE

function ParentCall {
### Script 03_run_parentcall.sh modified and simplified

    printf "Step 3 - Run ParentCall2.\n"

    #It is possible to map multiple families separately (different inversion rearrangement - marker order is not expected to be the same)
    #, I create here a loop with family name
    # note that lepmap allow several families to build the map together, they can all be put together in the pedigree file
    #cat $FAMILY | while read i
    #do
    #echo "family" $i

    #variables 
    #compulsory (file names)
    #GENO_FILE=$INPUT
    CALL_FILE="${OUTDIR}/03_parent_call_data/data_call_"$CORENAME".gz"
    # make the directory to be on the safe side
    mkdir -p "${OUTDIR}/03_parent_call_data"
    #touch $CALL_FILE
    #parameters
    #SEX="XLimit=2" #to call marker on sex chromosome in a XY system (use Zlimit=2 in a ZW system or nothing if we don't want sex-chromosome markers)

    #run the module
    #zcat $GENO_FILE | java -cp $LEPMAPDIR ParentCall2 data=$PEDIGREE posteriorFile=- removeNonInformative=1 | gzip > $CALL_FILE
    java -cp $LEPMAPDIR ParentCall2 data=$PEDIGREE vcfFile=$INPUT removeNonInformative=1 |gzip > $CALL_FILE
    #the parental call has also put the pedigree as header of the genotype data. To visualize:
    echo "these are the first 8 columns & 10 lines of the output file"
    zcat $CALL_FILE | cut -f 1-8 | head -n 10
    #done
}
ParentCall
printf "%-25s %s\n\n" "ParentCall2 output to:" $CALL_FILE >>$LOGFILE

function Filtering {
### Script 04_run_filtering apropriated and simplified 
    printf "\nStep 4 - Filtering\n"

    # this step will determine the markers actually used in the map.
    # the number id of each marker in the final map comes from the output of this stage

    # Read in variables 
    #parameters to adjust if needed
    read -p "Enter maximum missingness to be used in the filtering step - default [0.25]: " MISS
    MISS=${MISS:-0.25}
    echo "Maximum missingess for loci: $MISS"
    echo
    #MAF="MAFLimit=0.20" #unsure about how useful that parameter is for a single family : in one family min MAF should be 25% for all alleles
    read -p "Enter min MAF limit for use in the filtering step - default [0.05]: " MAF
    MAF=${MAF:-0.05}
    echo "MAF limit for loci: $MAF"
    echo
    #D="dataTolerance=0.001" #use a lower tolerance when running on several families
    read -p "Enter dataTolerance for filtering - default for single family crosses [0.0001]: " D
    D=${D:-0.0001}
    echo "dataTolerance=$D"
    echo

    #compulsory (file names)
    #IN_FILE="03_parent_call_data/data_call_"$MIN_COV"_"$i".gz" -- just use the CALL_FILE
    FILT_FILE="${OUTDIR}/04_filtering/data_f_"$CORENAME"_Miss-"$MISS".gz"
    CONTIG_FILE="${OUTDIR}/04_filtering/contig_"$CORENAME"_Miss-"$MISS".txt"
    POS_FILE="${OUTDIR}/04_filtering/pos_"$CORENAME"_Miss-"$MISS".txt"
    MARKERLIST_FILE="${OUTDIR}/04_filtering/data_f_"$CORENAME"_Miss-"$MISS".markerlist"
    #make sure the diectory exists
    mkdir -p "${OUTDIR}/04_filtering"

    #run the module
    zcat $CALL_FILE | java -cp $LEPMAPDIR Filtering2 data=- removeNonInformative=1 missingLimit=$MISS MAFLimit=$MAF dataTolerance=$D | gzip > $FILT_FILE
    #zcat ParentCall2_05.out.gz |java -cp /programs/LepMap3/bin/ Filtering2 data=- removeNonInformative=1 missingLimit=0.7 MAFLimit=0.05 dataTolerance=0.0001 |gzip >2_filtered_05.out.gz

    #extract the marker list
    echo "extracting the markerlist and formating it in R"
    zcat $FILT_FILE | awk '{print $1}' > $CONTIG_FILE
    zcat $FILT_FILE | awk '{print $2}' > $POS_FILE
    Rscript 00_scripts/Rscripts/4b_make_markerlizt.R "$CONTIG_FILE" "$POS_FILE" "$MARKERLIST_FILE" "$FAMILY"
    #done
    echo "Done filtering "
}
Filtering
{
    printf "%s\n" "Filtering parameters and output files"
    printf "%-25s %s\n" "Max missingness:" $MISS
    printf "%-25s %s\n" "MAF limit:" $MAF    
    printf "%-25s %s\n" "Data tolerance:" $D    
    printf "%-25s %s\n" "Filtering2 output:" $FILT_FILE 
    printf "%-25s %s\n" "Contig file:" $CONTIG_FILE
    printf "%-25s %s\n" "Position file:" $POS_FILE
    printf "%-25s %s\n" "Markerlist file:" $MARKERLIST_FILE
}>>$LOGFILE

###Step 5 Separate chromosomes - with some twists
#first make sure that the directory exists
mkdir -p "${OUTDIR}/05_map_chromosomes"
#Set default values before entering the loop
LOD=10 ; SIZEL=4

function SetLOD {
    #LOD=${LOD:-10}
    echo "LOD cutoff for SeparateChromosomes2 is currently set at [$LOD]"
    read -p "Enter new LOD cutoff value or <enter> to keep current setting : " NEWLOD
    LOD=${NEWLOD:-$LOD}
    echo "LOD cutoff set to $LOD"
}

function SetSIZEL {
    #SIZEL=${SIZEL:-4}
    echo "SIZEL cutoff for SeparateChromosomes2 is currently set at [$SIZEL]"
    read -p "Enter new SIZEL cutoff value or <enter> to keep current setting : " NEWSIZE
    SIZEL=${NEWSIZE:-$SIZEL}
    echo "SIZEL cutoff set to $SIZEL"
}

function OptimizeLod {
    printf "\nStep 5a: SeparateChromosomes - Optimize LOD \n"
    echo
    #echo "Now the important stuff begins - Pay attention !"
    read -p "Enter lower bound for LOD cutoff for SeperateChromosomes2 - default [5]: " LOD_LO
    LOD_LO=${LOD_LO:-5}
    echo $LOD_LO
    read -p "Enter upper bound for LOD cutoff for SeperateChromosomes2 - default [25]: " LOD_HI
    LOD_HI=${LOD_HI:-25}
    echo $LOD_HI
    echo
    echo "The minimum order of markers in order to call a linkage group (SIZEL) will be held constant while testing different LODs"
    SetSIZEL
    echo
    echo "Now we can run SeparateChromosomes2 repeatedly to find the optimal LOD cutoff."
    #variables 
    #compulsory (file names)
    LOD_LOG="${OUTDIR}/05_map_chromosomes/optimize_"$CORENAME"_Miss-"$MISS"_lod.log"

    printf "" >$LOD_LOG #clear old log the simple way
    printf "\n%s%s\n" "SIZEL fixed at " $SIZEL >$LOD_LOG.text
    
    for LODv in $(seq $LOD_LO 1 $LOD_HI)
    do
       echo -n "LOD $LODv "
       printf "LOD $LODv " >> $LOD_LOG
       #java -cp /usr/local/bin/lepmap3 SeparateChromosomes2 data=input_f.call lodLimit="${I}" sizeLimit=10 
       zcat $FILT_FILE | java -cp $LEPMAPDIR SeparateChromosomes2 data=- lodLimit=$LODv numThreads=$CPU sizeLimit=$SIZEL >/dev/null 2>>$LOD_LOG
    done
    echo
    #grep "Number of LGs" $LOD_LOG >$LOD_LOG.text
    grep -P 'LOD [0-9]{1,}|Number of LGs|^SIZEL' $LOD_LOG |sed -r 's/Loading file/:/'|sed '/:$/{N;s/\n/ /}' >>$LOD_LOG.text
    printf "Select the optimal LOD cutoff from the list\n"
    cat $LOD_LOG.text
}


function OptimizeSizel {
    printf "\nStep 5b: SeparateChromosomes - Optimize SIZEL \n"
    echo
    #echo "We have got the LOD now find the best SIZEL !"
    echo "Linkage groups with fewer than SIZEL loci will be discarded (but can later be rescued in the JoinSingle step)"
    read -p "Enter lower bound for SIZEL cutoff for SeperateChromosomes2 - default [2]: " SIZEL_LO
    SIZEL_LO=${SIZEL_LO:-2}
    echo $SIZEL_LO
    read -p "Enter lower bound for SIZEL cutoff for SeperateChromosomes2 - default [10]: " SIZEL_HI
    SIZEL_HI=${SIZEL_HI:-10}
    echo $SIZEL_HI
    SetLOD
    SIZEL_LOG="${OUTDIR}/05_map_chromosomes/optimize_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_sizel.log"
    
    printf "" >$SIZEL_LOG #clear old log the simple way
    printf "\n%s%s\n" "LOD fixed at " $LOD >$SIZEL_LOG.text

    for SIZELv in $(seq $SIZEL_LO 1 $SIZEL_HI)
    do
       echo -n "SIZEL $SIZELv "
       printf "SIZEL $SIZELv " >> $SIZEL_LOG
       #java -cp /usr/local/bin/lepmap3 SeparateChromosomes2 data=input_f.call lodLimit="${I}" sizeLimit=10 
       zcat $FILT_FILE | java -cp $LEPMAPDIR SeparateChromosomes2 data=- lodLimit=$LOD numThreads=$CPU sizeLimit=$SIZELv >/dev/null 2>>$SIZEL_LOG
    done
    echo
    #grep "Number of LGs" $SIZEL_LOG >$SIZEL_LOG.text
    grep -P 'SIZEL [0-9]{1,}|Number of LGs|^LOD' $SIZEL_LOG |sed -r 's/Loading file/:/'|sed '/:$/{N;s/\n/ /}' >>$SIZEL_LOG.text
    printf "Select the optimal SIZEL cutoff from the list\n"
    cat $SIZEL_LOG.text
}

while true; do
    printf "\n%s\n" "It is important to find the correct Lod cutoff and number of loci required to call a LG (Sizel)"
    printf "%s\n" "You can run iterations of either while keeping the other constant to select the best combination"
    printf "%s\n\n" "I will keep asking until you answer No"
    read -p "Do you want to optimize Lod or Sizel or Neither ? (l/s/n): " LSn
    case $LSn in
        [Ll]* ) OptimizeLod;;
        [Ss]* ) OptimizeSizel;;
        [Nn]* ) SetLOD; SetSIZEL; break;;
        * ) echo "Please answer L(od) S(izel) or N(o).";;
    esac
done

echo "Optimizations done"
{
    printf "\n%s\n" "Parameters for SeparateChromosomes2"
    printf "%-25s %s\n" "Selected LOD: " $LOD
    printf "%-25s %s\n" "Selected SIZEL: " $SIZEL
}>>$LOGFILE

function SeparateChromosomes {
    printf "\nStep 5c: Run SeparateChromosomes - Generate output files \n"
    echo "LOD set to $LOD and SIZEL set to $SIZEL"
    
    MAP_FILE="${OUTDIR}/05_map_chromosomes/map_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_Sizel-"$SIZEL".txt"
    REP_FILE="${OUTDIR}/05_map_chromosomes/map_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_Sizel-"$SIZEL".repartition"
    #other parameters (adjust if needed)
    ##### D="distortionLod=1 " #to remove markers with H-W distortion in a single family - !Z! this option is not being used as is.
    #RECOMB_1="maleTheta=0" # for species with non-recoùmbining male adjust recombination rate to 0
    #RECOMB_2="femaleTheta=0.03" # 
    #zcat $IN_FILE | java -cp bin/ SeparateChromosomes2 data=- lodLimit=$LOD numThreads=$CPU sizeLimit=$SIZEL $RECOMB_1 $RECOMB_2 > $OUT_FILE
    zcat $FILT_FILE | java -cp $LEPMAPDIR SeparateChromosomes2 data=- lodLimit=$LOD numThreads=$CPU sizeLimit=$SIZEL > $MAP_FILE

    #evaluate chromosome repartition 
    #typically if there is just one big chromosome -> raise LOD
    #if there are plenty of chromosome (more than expected -> lower LOD
    #if there are X big chromosome and plenty of chromosome with very few markers, raise SizeLimit and then rather use joinsingle with smalle LOD
    awk '{print $1}' $MAP_FILE | sort -n | uniq -c 
    awk '{print $1}' $MAP_FILE | sort -n | uniq -c > $REP_FILE
}
SeparateChromosomes

{
    printf "%-25s %s\n" "Map file: " $MAP_FILE
    printf "%-25s %s\n" "Repartition file: " $REP_FILE
}>>$LOGFILE

function JoinSingles {

    printf "\nStep 6: Run JoinSingles2All - Merge the singles into existing linkage groups \n"
    # The script this was based on had a "while" loop to iterate over multipe families - Could add this back later.

    mkdir -p "${OUTDIR}/06_join_singles"

    ### Prompt for different LOD and maybe other parameters
    echo "The current Lod cutoff is $LOD. You can use a lower Lod to merge singles into linkage groups"
    echo "this makes sense if there is a large number of markers not placed into LGs"
    read -p "Enter the Lod cutoff to use for JoinSingles2All - defaults to [$LOD]: " JSLOD
    JSLOD=${JSLOD:-$LOD}
    echo "Lod cutoff for JoinSingles2All set to $JSLOD"    
    #other parameters -- We don't need them now - Add later
    #D="distortionLod=1 " #to remove markers with H-W distortion in a single family
    #LODDIFF="lodDifference=1" #requires a difference of LOD # not understood: is it between LG?
    #RECOMB_1="maleTheta=0" # for species with non-recoùmbining male adjust recombination rate to 0
    #RECOMB_2="femaleTheta=0.03" # 
    #if there are many markers not placed into LG, this module can put them on the LG with a smaller LOD limit 
    cat <<ENDOFCOMMENT
This wrapper selects the most basic options for JoinSingles2All but there are a number of other options that can be useful:

Additional options:
         lodDifference=NUM  Required LOD difference [0.0]
         informativeMask=STR      Use only markers with informative father (1), mother(2), both parents(3) or neither parent(0) [0123]
         theta=NUM          Fixed recombination fraction [0.0 (Identicals)] [0.03 (All)]
         (fe)maleTheta=NUM  Fixed recombination fraction separately for both sex [theta]
                            set maleTheta=0 for species with non-recombining males
         betweenSameType=1  Only compute LOD scores between identically informative markers (applicable to single family data and Join...Identicals)
                            Also one has to specify 3 LOD limits (paternal, maternal and both) in this case
         lod3Mode=NUM       Controls how LOD scores are computed between double informative markers [1]
         maxDistance=NUM    Only calculate LOD scores between this many markers (JoinSingles2All only) [not set]
         distortionLod=1    Use segregation distortion aware LOD scores (JoinSingles2All only) [not set]
         mask=map_file2     Filter out markers not (in group) 1 in the map_file2

ENDOFCOMMENT
    read -p "Enter the extra parameters as shown above, separated by space  - default none: "
    EXTRAOPTIONS=$REPLY

    echo "Additional parameters to JoinSingles2All: $EXTRAOPTIONS"

    JS_MAP_FILE="${OUTDIR}/06_join_singles/js_map_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".txt"
    JS_REP_FILE="${OUTDIR}/06_join_singles/js_map_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".repartition"

    zcat $FILT_FILE | java -cp $LEPMAPDIR JoinSingles2All data=- lodLimit=$JSLOD numThreads=$CPU iterate=1 sizeLimit=$SIZEL map=$MAP_FILE $EXTRAOPTIONS > $JS_MAP_FILE
    TRY=$?
    #evaluate chromosome repartition 
    #typically if there is just one big chromosome -> raise LOD
    #if there are plenty of chromosome (more than expected -> lower LOD
    #if there are X big chromosome and plenty of chromosome with very few markers, raise SizeLimit and then rather use joinsingle with smaller LOD
    awk '{print $1}' $JS_MAP_FILE | sort -n | uniq -c 
    awk '{print $1}' $JS_MAP_FILE | sort -n | uniq -c > $JS_REP_FILE

    ### TODO add code to compare $REP_FILE and $JS_REP_FILE

    {
        printf "\n%s\n" "Parameters for JoinSingles2All"
        printf "%-25s %s\n" "LOD for JoinSingles2All: " $JSLOD
        printf "%-25s %s\n" "Additional parameters:" $EXTRAOPTIONS
        printf "%-25s %s\n" "Join Singles Map: " $JS_MAP_FILE
        printf "%-25s %s\n" "JS Repartition file: " $JS_REP_FILE
        if  [ ! $TRY -eq 0 ] ; then
            printf "JoinSingles2All - Command failed - something wrong with parameters\n"
        else
            printf "Success!"
        fi
    }>>$LOGFILE
}

TRY=1
JoinSingles
until [ $TRY -eq 0 ] ; do
    read -p "Do you want to try again? (y/n):" AGAIN
    case $AGAIN in
        [Yn]* ) JoinSingles;;
        [Nn]* ) echo "OK - skipping the step" ; break;;
        * ) echo "Please answer Yes or No.";;
    esac
done

function Order_markers {
    printf "\nStep 7: Run OrderMarkers2 and match_marker_map.R to ... \n"

    mkdir -p "${OUTDIR}/07_order_LG"
    #IN_FILE=$FILT_FILE
    #MAP_FILE="05_map_chr/map_"$MISS"_"$MIN_COV"_"$i"_LOD_"$LOD".txt" 
    #MAP_FILE="06_joinsingle/map.js_"$MISS"_"$MIN_COV"_"$i"_LOD_"$LOD".txt" # or use the map resulting from joinsingle 
    while true; do
        read -p "Use JoinSingles map file? (if not, the original map will be used) " yn
        case $yn in
            [Yy]* ) MAP=$JS_MAP_FILE; REP=$JS_REP_FILE; break;;
            [Nn]* ) MAP=$MAP_FILE; REP=$REP_FILE; break;;
            * ) echo "Please answer yes or no.";;
        esac
    done

    # NB_CHR=$[$(wc -l $REP | cut -d " " -f 1)-2]  ##No good on MacOS
    NB_CHR=$(tail -1 $REP | awk '{print $2}')  #number of linkage groups from step 5
    printf "\nNumber of Chromosomes: $NB_CHR\n"

    #other parameters (adjust if needed)
    read -p "How many iterations shall we use for merging? - default [5]: " ITE
    ITE=${ITE:-5}
    echo "numMergeIterations=$ITE"
    read -p "How many Refine steps shall we use? - default [2]: " REFINE_STEPS
    REFINE_STEPS=${REFINE_STEPS:-2}
    echo "Refine steps=$REFINE_STEPS"

    #RECOMB_1="recombination1=0" # for species with non-recoùmbining male adjust recombination rate to 0
    #PHASE="outputPhasedData=1" #if we want phased data as output. may be useful for QTL?

    #run the module once for intitial order
    for j in $(seq $NB_CHR)
    do
        printf "\n* Assessing marker order for LG$j \n\n"
        OUT_FILE="${OUTDIR}/07_order_LG/order_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".LG"$j".1.txt"
        echo "INPUT: $FILT_FILE"

        zcat $FILT_FILE | java -cp $LEPMAPDIR OrderMarkers2 map=$MAP numThreads=$CPU data=- numMergeIterations=$ITE chromosome=$j $RECOMB_1 $PHASE > $OUT_FILE

        for k in $(seq $REFINE_STEPS)
        do
            IT=$[$k + 1]                           #Everything is global in bash so we can use this later in LMplots sub
            printf "\n** Assessing marker order for LG$j refining step $IT \n\n"
            OUT_FILE="${OUTDIR}/07_order_LG/order_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".LG"$j"."$IT".txt"
            ORDER_FILE="${OUTDIR}/07_order_LG/order_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".LG"$j"."$k".txt"

            zcat $FILT_FILE | java -cp $LEPMAPDIR OrderMarkers2 evaluateOrder=$ORDER_FILE numThreads=$CPU data=- numMergeIterations=$ITE chromosome=$j $RECOMB_1 $PHASE > $OUT_FILE
        done
    done
}
Order_markers

{
    printf "\n%s\n" "OrderMarkers2 parameters"
    printf "%-25s %s\n" "Number of LGs: " $NB_CHR
    case $yn in
        [Yy]* ) echo "JoinSingles Map used";;
        [Nn]* ) echo "Original Map file used";;
    esac
    printf "%-25s %s\n" "Output written to: " "${OUTDIR}/07_order_LG/order_*"
}>>$LOGFILE


function match_markers_and_plot {
    printf "\nStep 8: Analyzing maps and making some plots\n"

    mkdir -p "${OUTDIR}/08_analyze_maps/01_maps"
    mkdir -p "${OUTDIR}/08_analyze_maps/02_plots"
    ORDERFILES="${OUTDIR}/07_order_LG/order_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".LG"
    OUTFILEBASENAME=$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL
    OUTDIRBASE="${OUTDIR}/08_analyze_maps/"

    #This Rscript uses the marker list from 04_filtering to assign a contig & pos to the marker id given by lepmap,
    #it makes the map more easily explorable, and correspondance with the genome used to align the markers
    #last step also outputs a plot with the female position which allow looking at the repartition of markers along each LG

    Rscript 00_scripts/Rscripts/8_match_marker_map_fx.R $MARKERLIST_FILE $NB_CHR $REFINE_STEPS $ORDERFILES $OUTFILEBASENAME $OUTDIRBASE

}
match_markers_and_plot

function Make_LMPlots {

    printf "\nMaking some nice plots to view with xdot.py - placing them in: 08_analyze_maps/03_LMPlots\n"
    INPUT="${OUTDIR}/07_order_LG/order_"$CORENAME"_Miss-"$MISS"_Lod-"$LOD"_JSLod-"$JSLOD"_Sizel-"$SIZEL".LG"
    #i.e. the orderfiles from order markers

    mkdir -p "${OUTDIR}/08_analyze_maps/03_LMPlots"

    for LG in $(seq $NB_CHR)
    do
        #echo "Inputfile: 07_order_LG/"$INPUTFILENAME".LG"$LG".3.txt"
        #echo "Outputfile: 09_LMplots/"$INPUTFILENAME".LG"$LG".3.dot"
        java -cp $LEPMAPDIR LMPlot $INPUT""$LG".$IT.txt" > "${OUTDIR}/08_analyze_maps/03_LMPlots/"${OUTFILEBASENAME}".LG"${LG}".$IT.dot"  #IT is the highest iteration 
    done
}
Make_LMPlots

function Fmt_Mpcmp_and_Chrmnmr {

    function Get_genome_and_run {
        read -e -p "Enter the path to the indexed genome.fasta file (can be bgzipped): " GENOME

        printf "command:  bash ./00_scripts/01_Extract_loci-format_for_mapcomp_and_chromonomer.sh $GENOME $OUTDIR\n"
        bash ./00_scripts/01_Extract_loci-format_for_mapcomp_and_chromonomer.sh $GENOME $OUTDIR
        
        if [[ $? = 0 ]] ; then
            DONE=1
            printf "This seems to have worked \n"
            printf "Sucessfully formatted files for mapcomp and chromonomer.\n">>$LOGFILE
        else
            printf "Something didn't work. Run script again or run it later on it's own.\n\n"
            printf "Tried to format for mapcomp and chromonomer but failed.\n">>$LOGFILE
            AG="again "
        fi
    }

    printf "\nOLD STYLE: Optional formatting for Mapcomp and Chromonomer input (up to a point).\n\n"
    printf "The outputfiles can be useful even if you do not plan on using either.  You will have to provide\n"
    printf "the path to the indexed fasta file which was used for alignment when SNPs were called\n"
    printf "If snps were called \"denovo\" unsing Stacks you should skip this step.  It won't get you anywhere\n"
    printf "\nOutput will be written subdirectories of 08_analyze_maps/ (old style)\n"

    while true; do
        read -p "Do you want to run the script ${AG}? (y/n): " RUNIT
        case $RUNIT in
            [Yy]* ) echo "Yes!"; Get_genome_and_run ; if [[ $DONE = 1 ]]; then break ; fi ;;
            [Nn]* ) echo "OK - skipping the step" ; printf "Skipped making Mapcomp and Chromonomer files\n" ; break;;
            * ) echo "Please answer Yes or No.";;
        esac
    done
}

function samtools_genome_extractor {

    MAPDIR="${OUTDIR}/08_analyze_maps/01_maps/"
    SEQOUTDIR="${OUTDIR}/09_extracted_genome_seqs/"
    mkdir -p $SEQOUTDIR
    
    printf "\nReading inputfiles from:\t$MAPDIR\n"
    printf "Writing output files to:\t$SEQOUTDIR\n"
    
    printf "\nSamtools should be installed or this will not work.\n\n"
    printf "You must provide the path to the indexed fasta file which was used for alignment to call SNPs\n"
    while true; do
        read -e -p "Enter the path to the indexed genome.fasta file (can be bgzipped): " GENOME
        if [ -f "$GENOME" ]; then
            printf "\n$GENOME file found.\n"
            break
        else 
            echo "file $GENOME does not exist. Please try again!"
            read -p "Do you want to try again or quit ? (y/q): " RUNIT
            case $RUNIT in
                [Yy]* ) echo "Yes!";;
                [Qq]* ) echo "OK - giving up" ; printf "Skipped making Mapcomp and Chromonomer files\n" ; exit 1;;
                * ) echo "Please answer Yes or Quit.";;
            esac
        fi
    done
    
    ls $MAPDIR | while read i
    do 
         echo "Infile: $i"
         printf "command: perl ../LepMapWrapper/00_scripts/01b_genome_extractor.pl -m ${MAPDIR}${i} -g $GENOME -o $SEQOUTDIR\n"
         printf "\n"
         (perl ./00_scripts/01b_genome_extractor.pl -m ${MAPDIR}${i} -g $GENOME -o $SEQOUTDIR </dev/tty)
         echo "returned from perl"
    done
}

function catalog_extractor {

    MAPDIR="${OUTDIR}/08_analyze_maps/01_maps/"
    SEQOUTDIR="${OUTDIR}/09_extracted_genome_seqs/"
    mkdir -p $SEQOUTDIR
    
    while true; do
        read -e -p "Enter the path to the stacks catalog.fasta.gz file : " CATALOG
        if [ -f "$CATALOG" ]; then
            printf "\n$CATALOG file found.\n"
            denovo_pattern="^>[0-9]*\sNS=[0-9]*\scontig="
            refmap_pattern="^>[0-9]*\spos=[[:alnum:]]*.*NS=[0-9]*"
            firstline=$(zcat $CATALOG | head -n 1)
            printf "First line of catalog file: $firstline  => "
            if [[ $firstline =~ $denovo_pattern ]]; then
                catalogtype="denovo"
            elif [[ $firstline =~ $refmap_pattern ]]; then
                catalogtype="reference_aligned"
            else
                catalogtype="undefined"
            fi
    
            printf "Catalog_type detected as: $catalogtype\n\n"
            if [[ $catalogtype =~ "undefined" ]]; then echo "This does not look good - quitting now\n" ; exit 1; fi
            break
        else 
            echo "file $CATALOG does not exist. Please try again!"
            read -p "Do you want to try again or quit ? (y/q): " RUNIT
            case $RUNIT in
                [Yy]* ) echo "Yes!";;
                [Qq]* ) echo "OK - giving up" ; printf "Skipped making Mapcomp and Chromonomer files\n" ; exit 1;;
                * ) echo "Please answer Yes or Quit.";;
            esac
        fi
    done
    
    NINETYFIVE=$((NUM_SAMPLES*95/100))
    printf "\nDiscarding low coverage stacks from the catalog will seriously speed up things.\n"
    printf "A sample number close to 95% of the number of individuals in family is recommended."
    read -p "Enter minimal samplenumber threshold - suggested [$NINETYFIVE]: " DISCARD
    DISCARD=${DISCARD:-$NINETYFIVE}
    
    ls $MAPDIR | while read i
    do 
         echo "Infile: $i"
         if [[ $catalogtype =~ "reference_aligned" ]]; then
            printf "command:  perl ./00_scripts/02_catalog_extractor.pl -d $DISCARD -m ${MAPDIR}${i} -c $CATALOG -o $SEQOUTDIR \n"
            (perl ./00_scripts/02_catalog_extractor.pl -m ${MAPDIR}${i} -c $CATALOG -o $SEQOUTDIR -d $DISCARD </dev/tty)
         elif [[ $catalogtype =~ "denovo" ]]; then
            printf "command:  perl ./00_scripts/02b_denovo_catalog_extractor.pl -d $DISCARD -m ${MAPDIR}${i} -c $CATALOG -o $SEQOUTDIR \n"
            (perl ./00_scripts/02b_denovo_catalog_extractor.pl -m ${MAPDIR}${i} -c $CATALOG -o $SEQOUTDIR -d $DISCARD </dev/tty)            
         fi
         echo "catalog extracted"
    done
}

function Ask_old_or_new {
    while true; do
        printf "\n%s\n\n" "OK, there are two ways to go now, using a genome:"
        printf "%s\n" "  1) Use python and R scripts from 'clairemerot' - (old style)"
        printf "%s\n\n" "  2) Use samtools to extract sequences - (faster)"
        read -p "Old style or New or Quit ? (o/n/q) : " OldOrNew
        case $OldOrNew in
            [Oo]* ) Fmt_Mpcmp_and_Chrmnmr ; break;;
            [Nn]* ) samtools_genome_extractor ; break;;
            [Qq]* ) break;;
            * ) echo "Please answer G(enome) S(tacks) or Q(uit).";;
        esac
    done
}

while true; do
    printf "\n\n%s\n\n" "Step 9: Extracting sequences flanking the markers to use for further analysis"
    printf "%s\n" "We have different ways of searching for sequences and preparing output files that"
    printf "%s\n" "can be used for input for tools such as Mapcomp and Chromonomer (after some tweaking)."
    printf "%s\n" "Depending on how you got your markers you can extract sequences from the genome that you"
    printf "%s\n" "used to call your SNPs or from a stacks catalog.gz file."
    printf "\n%s\n" "You have to make some choices now:"
    read -p "Use Genome or Stacks catalog or just call it Quits ? (g/s/q) : " GSQ
    case $GSQ in
        [Gg]* ) Ask_old_or_new ; break;;
        [Ss]* ) catalog_extractor ; break;;
        [Qq]* ) break;;
        * ) echo "Please answer G(enome) S(tacks) or Q(uit).";;
    esac
done

while true; do
    printf "\n%s\n\n" "Do you want to generate more output files?"
    printf "%s\n\n" "Same options as before. Genome (old / new scripts) or Stacks catalog"
    read -p "Extract regions from genome using 'New' or 'Old' or 'Stacks catalog' or Quit ? (n/o/s/q) : " GSQagain
    case $GSQagain in
        [Nn]* ) samtools_genome_extractor ;;
        [Oo]* ) Fmt_Mpcmp_and_Chrmnmr ;;
        [Ss]* ) catalog_extractor ;;
        [Qq]* ) break;;
        * ) echo "Please answer G(enome) S(tacks) or Q(uit).";;
    esac
done
    

echo "End of the pipe - Inspect the output files and log to see whether it ran successfully!"

#TODO Add more logging towards the "end of the line"
