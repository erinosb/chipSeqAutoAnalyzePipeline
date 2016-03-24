#! /bin/sh/
##################################################################################################

#####################  HEADER  ###################################################################
#  chipSeqAnalyzeStep1.sh
#  Copyright (c) 2016, Erin Osborne Nishimura
#
#PROGRAM
#   chipSeqAnalyzeStep1.sh - To automate the quality control & alignment of multiplexed ChIP-seq data
#
#
#USAGE
#   Split and Align Mode:
#       bash chipSeqAnalyzeStep1.sh [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
#   OR
#   Split Only Mode:
#       bash chipSeqAnalyzeStep1.sh --splitOnly [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
#   OR
#   Align Only Mode:
#       bash chipSeqAnalyzeStep1.sh --alignOnly [options] <input1.fastq> input2.fastq
#
#
#MODES
#   There are three modes for use... --splitNAlign mode (default), --splitOnly and --alignOnly mode
#   --splitNAlign (default)      This mode uses as input a single homebrew multiplexed .txt file and
#                               a single barcode index file. The pipeline then...
#                                   0) splits the multiplexed file into single sample files by index
#                                   1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
#                                   2) perform Tagdust using a solexa library of adapter and primer sequences.
#                                   3) perform Fastqc to make a report of quality
#                                   4) bowtie alignment
#                                   5) compress .sam --> .bam using samtools view
#                                   6) Sort .bam into _sorted.bam using samtools sort
#                                   7) Convert _sorted.bam into .bed using bedtools
#                                   8) Convert .bed into .bw file using bedToBw.sh dependency script
#    
#   --splitOnly                 This mode uses as input a single homebrew multiplexed .txt file and
#                               a single barcode index file. The pipeline then...
#                                   0) splits the multiplexed file into single sample files by index
#
#   --alignOnly                  This mode uses as input a list of .txt or .fastq files and analyzes them...
#                               a single barcode index file. The pipeline then...
#                                   1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
#                                   2) perform Tagdust using a solexa library of adapter and primer sequences.
#                                   3) perform Fastqc to make a report of quality
#                                   4) bowtie alignment
#                                   5) compress .sam --> .bam using samtools view
#                                   6) Sort .bam into _sorted.bam using samtools sort
#                                   7) Convert _sorted.bam into .bed using bedtools
#                                   8) Convert .bed into .bw file using bedToBw.sh dependency script
#
#
#ARGUMENTS
#   <inputFile.txt>           This is the file from the sequencing facility. It is a fastq file containing Illumina sequencing reads generated from multiplexed samples
#   <barcodeIndexFile.txt>    This is a file listing the barcodes that are at the 5' end of each sequence string.  It should be in the format:
#
#                   #barcode file for the multiplexed sequences generated 11/11/12
#                   EO33	AGATGGT
#                   EO34	CACGTCG
#                   EO35	GATCTTG
#                   EO36	TCAGGAC
#                   EO37	ACAGTTG
#   <input1.fastq>          This is an input .fastq file 
#
#
#OPTIONS
#
#   --splitOnly                      Runs in splitOnly mode. Suppresses analysis.
#   --alignOnly                      Runs in alignOnly mode. Suppresses splitting sequencefiles
#   --qualityOff                     Suppresses quality score reports
#   --cleanOff                       Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
#                                      It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
#   --trimOff                        Suppresses trimming six basepairs of multiplexing indices
#
#   --extension <n>                  Specify the length bp to extend reads in the .wig file. Default = 100
#   -p <n>                           Runs bowtie in parallel mode. Values accepted are 1 - 8. Suggest 2 - 4. Make sure to match this number with bsub -n <n>. Default = 1
#   --chrlength <file.txt>           A file specifying how long each chromsome is. The bedToBw.sh program needs this information for read extension. Default =
#                                       /proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt
#   --bowtiepath </path/>            The path location of where the .bwa files are contained. Required by bowtie. Default =
#                                       /proj/dllab/Erin/ce10/from_ucsc/seq/prev_versions_bowtie/genome_bwa/ce10
#   --primers <file.txt>             A .fastq file that contains a list of all the primer and adpater sequences used in library prep and sequencing. Tagdust
#                                      remove these sequences as a quality control step. Default =
#                                       /proj/dllab/Erin/sequences/solexa-library-seqs.fasta
#   
#
#DEPENDENCIES
#   Requires tagdust, fastqc, fastx-toolkit, bowtie1, bedtools, samtools v1.0 or greater, bedGraphToBigWig
#   Developed with versions: TagDust 1.12; FastQC v0.11.3; bowtie/1.1.2; bedtools/2.25.0; samtools/1.3
#
#
#BUGS/FUTURE EXPANSION
   #-- check whether certain modules have been loaded, exit with an error if they are missing.
   #-- Auto load all required modules
   #-- toggle between bowtie or bowtie2
   #-- add an option for an output directory --outdir
   #-- add an option for a 'nickname' for the logfiles --logfile
   
# Updates
    # 160323 --> Updated new samtools -o option
    # 160323 --> added a chunkmem option to bowtie to avoid warnings
    # 160323 --> Fixed bgs with AlignOnly method
    # 150323 --> Fixed bugs with summary
#
#
#
##################################################################################################

#####################   SET VARIABLES   ##########################################################
solexa_primer_adapter="/proj/dllab/Erin/sequences/solexa-library-seqs.fasta"    #tagdust needs a .fasta file that contains a list of all the solexa primer and adapter sequences. Set this
                                                                                  #variable to a path pointing to that file.
#bowtie2path="/proj/dllab/Erin/ce10/from_ucsc/seq/genome_bt2/ce10"               #bowtie2 know where the bowtie2 index files are located. Set this varaible to the path and root
                                                                                  #name of those index files.
                                                                                  #Also, the genome sequence (a .fa file) also needs to be in that same directory.
bowtie1path="/proj/dllab/Erin/ce10/from_ucsc/seq/prev_versions_bowtie/genome_bwa/ce10"
chromlength="/proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt"               #bedToBw.sh needs to know how long each chromosome is
##################################################################################################

#####################  USAGE   ###################################################################

usage="  chipSeqAnalyzeStep1.sh
  Copyright (c) 2015, Erin Osborne Nishimura

PROGRAM
   chipSeqAnalyzeStep1.sh - To automate the quality control & alignment of multiplexed ChIP-seq data


USAGE
   Split and Align Mode:
       bash chipSeqAnalyzeStep1.sh [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
   OR
   Split Only Mode:
       bash chipSeqAnalyzeStep1.sh --splitOnly [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
   OR
   Align Only Mode:
       bash chipSeqAnalyzeStep1.sh --alignOnly [options] <input1.fastq> input2.fastq


MODES
   There are three modes for use... --splitNAlign mode (default), --splitOnly and --alignOnly mode
   --splitNAlign (default)      This mode uses as input a single homebrew multiplexed .txt file and
                               a single barcode index file. The pipeline then...
                                   0) splits the multiplexed file into single sample files by index
                                   1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
                                   2) perform Tagdust using a solexa library of adapter and primer sequences.
                                   3) perform Fastqc to make a report of quality
                                   4) bowtie alignment
                                   5) compress .sam --> .bam using samtools view
                                   6) Sort .bam into _sorted.bam using samtools sort
                                   7) Convert _sorted.bam into .bed using bedtools
                                   8) Convert .bed into .bw file using bedToBw.sh dependency script
    
   --splitOnly                 This mode uses as input a single homebrew multiplexed .txt file and
                               a single barcode index file. The pipeline then...
                                   0) splits the multiplexed file into single sample files by index

   --alignOnly                  This mode uses as input a list of .txt or .fastq files and analyzes them...
                               a single barcode index file. The pipeline then...
                                   1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
                                   2) perform Tagdust using a solexa library of adapter and primer sequences.
                                   3) perform Fastqc to make a report of quality
                                   4) bowtie alignment
                                   5) compress .sam --> .bam using samtools view
                                   6) Sort .bam into _sorted.bam using samtools sort
                                   7) Convert _sorted.bam into .bed using bedtools
                                   8) Convert .bed into .bw file using bedToBw.sh dependency script


ARGUMENTS
   <inputFile.txt>           This is the file from the sequencing facility. It is a fastq file containing Illumina sequencing reads generated from multiplexed samples
   <barcodeIndexFile.txt>    This is a file listing the barcodes that are at the 5' end of each sequence string.  It should be in the format:

                   #barcode file for the multiplexed sequences generated 11/11/12
                   EO33	AGATGGT
                   EO34	CACGTCG
                   EO35	GATCTTG
                   EO36	TCAGGAC
                   EO37	ACAGTTG
   <input1.fastq>          This is an input .fastq file 


OPTIONS

   --splitOnly                      Runs in splitOnly mode. Suppresses analysis.
   --alignOnly                      Runs in alignOnly mode. Suppresses splitting sequencefiles
   --qualityOff                     Suppresses quality score reports
   --cleanOff                       Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
                                      It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
   --trimOff                        Suppresses trimming six basepairs of multiplexing indices

   --extension <n>                  Specify the length bp to extend reads in the .wig file. Default = 100
   -p <n>                           Runs bowtie in parallel mode. Values accepted are 1 - 8. Suggest 2 - 4. Make sure to match this number with bsub -n <n>. Default = 1
   --chrlength <file.txt>           A file specifying how long each chromsome is. The bedToBw.sh program needs this information for read extension. Default =
                                       /proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt
   --bowtiepath </path/>            The path location of where the .bwa files are contained. Required by bowtie. Default =
                                       /proj/dllab/Erin/ce10/from_ucsc/seq/prev_versions_bowtie/genome_bwa/ce10
   --primers <file.txt>             A .fastq file that contains a list of all the primer and adpater sequences used in library prep and sequencing. Tagdust
                                      remove these sequences as a quality control step. Default =
                                       /proj/dllab/Erin/sequences/solexa-library-seqs.fasta
   

DEPENDENCIES
   Requires tagdust, fastqc, fastx-toolkit, bowtie1, bedtools, samtools v1.0 or greater, bedGraphToBigWig
   Developed with versions: TagDust 1.12; FastQC v0.11.3; bowtie/1.1.2; bedtools/2.25.0; samtools/1.3"

##################################################################################################

#####################
#PRE-PROCESSING: Start Log files, load modules, Check errors, get filenames, set options
#####################

#Start logfiles
DATE=$(date +"%Y-%m-%d_%H%M")
dated_log=${DATE}.log
commands_log=${DATE}_commands.log

printf "\n\n"  | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log
printf $(date +"%Y-%m-%d_%H:%M")"\tRUNNING\n" | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log


printf $(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log

echo -e "INITIATED autoAnalyzeChipseq_v9.sh using command:\n\t$0 $*\n" | tee -a $dated_log $commands_log

printf "This pipeline requires the following modules: samtools, bedtools, bowtie, fastqc"

printf "\nThis pipeline will run with the following modules:\n" | tee -a $dated_log $commands_log ##Doesn't work
source /nas02/apps/Modules/default/init/bash
module list  2>&1 | tee -a $dated_log $commands_log
printf "\n" | tee -a $dated_log $commands_log


#Set options
alignonly="notcalled"
qualityoff="notcalled"
splitonly="notcalled"
trimoff="notcalled"
multi="undefined"
bar="undefined"
parallel=1
extension=100
cleanoff="notcalled"
list=$@

#Check for options and input:
if [ -z "$1" ]
  then
    echo "ERROR: No options or input files supplied:" | tee -a $dated_log $commands_log
    echo "$usage"  | tee -a $dated_log $commands_log
    exit
fi

#Get options
for n in $@
do
    case $n in
        --splitNAlign) shift;
        ;;
        --alignOnly) shift;
        alignonly="called"
        ;;
        --splitOnly) shift;
        splitonly="called"
        ;;
        --trimOff) shift;
        trimoff="called"
        ;;
        --qualityOff) shift;
        qualityoff="called"
        ;;
        --extension) shift;
        extension=${1:--}
        shift
        ;;
        --multi) shift;
        multi="$1"
        shift
        ;;
        --bar) shift;
        bar="$1"
        shift
        ;;
        -p) shift;
        parallel=${1:--}
        shift
        ;;
        --cleanOff) shift;
        cleanoff="called"
        ;;
        --primers) shift;
        solexa_primer_adapter=${1:--}
        shift
        ;;
        --bowtiepath) shift;
        bowtie1path=${1:--}
        shift
        ;;
        --chrlength) shift;
        chromlength=${1:--}
        shift
        ;;
    esac
done

#Check for errors. If split is to be performed, there should be a multi file and a bar file:

#Report that splitonly mode was called. Check for errors if splitonly mode is called.
if [[ $splitonly == "called" ]]
then
    printf "\nRUNNING IN SPLITONLY MODE\n"  | tee -a $dated_log $commands_log
    printf "\nFILE TO SPLIT IS: $multi" | tee -a $dated_log $commands_log
    printf "\nBARCODE FILE IS: $bar" | tee -a $dated_log $commands_log
    
    if [[ $multi == "undefined" ]]
    then
        printf "\n" | tee -a $dated_log $commands_log
        echo "ERROR: No inputFile supplied:" | tee -a $dated_log $commands_log
        echo "$usage"  | tee -a $dated_log $commands_log
        exit
    fi
    
    if [[ $bar == "undefined" ]]
    then
        printf "\n" | tee -a $dated_log $commands_log
        echo "ERROR: No barcodeIndexFile supplied:" | tee -a $dated_log $commands_log
        echo "$usage"  | tee -a $dated_log $commands_log
        exit
    fi
fi

#Report that splitnalign mode was called. Check for errors if splitnalign mode is called.
if [[ $alignonly == "notcalled" && $splitonly == "notcalled" ]]
then
    printf "\nRUNNING IN SPLIT-N-ALIGN MODE\n"  | tee -a $dated_log $commands_log
    printf "\nFILE TO SPLIT IS: $multi" | tee -a $dated_log $commands_log
    printf "\nBARCODE FILE IS: $bar" | tee -a $dated_log $commands_log
    
    if [[ $multi = "undefined" ]]
    then
        printf "\n" | tee -a $dated_log $commands_log
        echo "ERROR: No inputFile supplied:" | tee -a $dated_log $commands_log
        echo "$usage"  | tee -a $dated_log $commands_log
        exit
    fi
    
    if [[ $bar = "undefined" ]]
    then
        printf "\n" | tee -a $dated_log $commands_log
        echo "ERROR: No barcodeIndexFile supplied:" | tee -a $dated_log $commands_log
        echo "$usage"  | tee -a $dated_log $commands_log
        exit
    fi
 
fi

#Report that alignonly mode was called. Check for errors if alignonly mode is called.
if [[ $alignonly == "called" ]]
then
    printf "\nRUNNING IN ALIGNONLY MODE\n"  | tee -a $dated_log $commands_log
    printf "\nFILES TO ALIGN ARE: $*" | tee -a $dated_log $commands_log
    
    if [ -z "$1" ]
    then
        printf "\n" | tee -a $dated_log $commands_log
        echo "ERROR: No options or input files supplied:" | tee -a $dated_log $commands_log
        echo "$usage"  | tee -a $dated_log $commands_log
        exit
    fi
fi

######################
#If splitting needs to happen, un-multiplex using fastx barcode splitter [alignOnly == "notcalled"]
######################

if [[ $alignonly == "notcalled" ]]
then
    printf "\n\n"  | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log
    printf $(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "SPLITTING" | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log

    #printf "\n\tSplitting $multi into multiple files based on barcodes in $bar\n"  | tee -a $dated_log $commands_log
    #more $bar | tee -a $dated_log $commands_log
    
    
    ######################
    #Split the multiplexed file into multiple files based on barcoded indexes
    ######################
    printf "\n\t\t"
    if [[ ${multi} == *.gz ]]
    then
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\tFASTX Toolkit: Splitting multiplexed fastq file based on barcode indices.\n\t"
        echo "zcat $multi | fastx_barcode_splitter.pl --bcfile $bar --prefix "" --suffix ".fastq" --bol" | tee -a $dated_log $commands_log
        zcat $multi | fastx_barcode_splitter.pl --bcfile $bar --prefix "" --suffix ".fastq" --bol | tee -a $dated_log
    else
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\tFASTX Toolkit: Splitting multiplexed fastq file based on barcode indices.\n\t"
        echo "cat $multi | fastx_barcode_splitter.pl --bcfile $bar --prefix "" --suffix ".fastq" --bol" | tee -a $dated_log $commands_log
        cat $multi | fastx_barcode_splitter.pl --bcfile $bar --prefix "" --suffix ".fastq" --bol | tee -a $dated_log        
    fi
    
fi




##############################
#If splitOnly is called, exit
##############################

if [[ $splitonly == "called" ]]
then
    exit

fi
    
    
    
#############################
#Generate a list of files that need to be aligned and analyzed in the next code block.
#############################

#list is an array that contains the prefixes of all the .fastq files to analyze.
list=()
fileextension=
#path=

#If alignonly is called, the .fastq files are still in $@
if [[ $alignonly == "called" ]]
then
    input_files=$*
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\t" | tee -a $dated_log $commands_log
    printf "Will Chip-seq analyze the following .fastq/.txt files:\n" | tee -a $dated_log $commands_log
    for i in $input_files
        do
            echo ${i} | tee -a $dated_log  $commands_log
        done
    
    #Remove the .fastq or .txt postfix from the name
    for i in $input_files
        do
            #Parse input file:
            path=$(dirname "${i}")
            #echo "path is $path"
    
            myfile=$(basename $i)
            #echo "myfile is $myfile"  | tee -a $dated_log  $commands_log
            root="${myfile%%.*}"
            #echo "my root is $root" | tee -a $dated_log  $commands_log
            list+=("${root}")
            #echo "inside loop list is ${list[@]}" | tee -a $dated_log  $commands_log
        done
        
        
    ##Get file extension
    #firstfile="${input_files}"
    #fileextension="${firstfile##*.}"

    
fi

#In splitNAlign mode, the names of the files to analyze should be in the --bar barcode file
if [[ $alignonly == "notcalled" && $splitonly == "notcalled" ]]
then
    list=($(grep "\#" -v ${bar} | awk '{print $1}'))
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\t" | tee -a $dated_log $commands_log
    printf "Will Chip-seq analyze the following .fastq/.txt files:\n" | tee -a $dated_log $commands_log
    for i in "${list[@]}"
        do
            echo ${i}".fastq" | tee -a $dated_log  $commands_log
        done
    fileextension="fastq"
fi


#############################
#Quality Processing, Alignment, Output
#############################

#For each sample,
#    1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
#    2) perform Tagdust using a solexa library of adapter and primer sequences.
#    3) perform Fastqc to make a report of quality
#    4) bowtie alignment
#    5) compress .sam --> .bam using samtools view
#    6) Sort .bam into _sorted.bam using samtools sort
#    7) Convert _sorted.bam into .bedgraph (.wig file) using zinba
  
input_files=()

if [[ $alignonly == "called" ]]
then
    input_files=$*
    #echo "alignonly called; input file array is @input_files"
else
    input_files=( "${list[@]/%/.fastq}" )
    #echo "alignonly not called; input file array is ${input_files[@]}"
fi


for j in ${input_files[@]}
do
    
    #echo "within the quality processing, alignment, output, j is $j"
    
    #Parse input file:
    path=$(dirname "${j}")
    
    myfile=$(basename $j)
    
    #Remove the suffix
    i="${myfile%%.*}"
    
    #echo "Processing: $i" | tee -a $dated_log  $commands_log   
    
    
    
    opd=${i}_opd
    opdpath="${i}_opd/"
    trimfile="1_${i}_trim.fastq"
    cleanfile="2_${i}_clean.fastq"
    samfile="4_${i}_output.sam"
    unaligned="4_${i}_unaligned.fastq"
    metrics="4_${i}_bowtiemetrics.txt"
    bamfile="5_${i}.bam"
    bam_sorted="6_${i}_sorted"
    bam_sort_file="6_${i}_sorted.bam"
    bed_file="7_${i}.bed"
    bw_file="8_${i}x${extension}n.bw"
    
    
    
    
    printf "\n\n"  | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log
    printf $(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "PROCESSING $i:" | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log
    
    #Make Output Directory
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "mkdir: Make output directory.\n\tmkdir $opd" | tee -a $dated_log $commands_log
    mkdir $opd 2>&1 | tee -a $dated_log

        #echo "trimOff setting is $trimoff"
        #echo "alignOnly setting is $alignonly"
        
        
    ## OK, this is complicated. If SplitAlign mode is called, then the reads were split and we need to trim them.
    ## If align only is called, the reads aren't trimmed.
    ## I guess in certain cases, one could trim and align only read (for example if the sequencer split the reads but did not trim them)
    ## Loop throug these options and parse out the different strategies. This is complicated by the fact that tagdust is either going to be working on a _trim file or on a .fastq file
    
    # splitNAlign mode
    if [ $trimoff == "notcalled" ] && [ $alignonly == "notcalled" ]
    then
        echo "Running in splitNAlign mode"
        #1) fastx_trim:
        cmd1="fastx_trimmer -f 9 -Q 33 -i "$i"."$fileextension" -o "$opdpath$trimfile
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "fastx_trimmer: Trimming barcode indexes from ${i}.${fileextension} using command:\n\t$cmd1" | tee -a $dated_log $commands_log
        fastx_clipper -h | grep 'FASTX' - 2>&1 | tee -a $dated_log $commands_log
        $cmd1  2>&1 | tee -a $dated_log
        
        #2) tagdust
        cmd2="tagdust -q -f 0.001 -s -a "$opdpath"2_"$i"_artifact.fastq -o "$opdpath$cleanfile" "$solexa_primer_adapter" "$opdpath$trimfile
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "Tagdust: Removing adapter and primer sequences from $trimfile to make $cleanfile using command:\n\t$cmd2" | tee -a $dated_log $commands_log
        $cmd2  2>&1 | tee -a $dated_log
    
    # The weird alignOnly mode with trimming
    elif [ $trimoff == "notcalled" ] && [ $alignonly == "called" ]
    then
        echo "Running in alignOnly mode but with trimming. Make sure you want to do this."
        
        #1) fastx_trim:
        cmd1="fastx_trimmer -f 9 -Q 33 -i "$j"."$fileextension" -o "$opdpath$trimfile
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "fastx_trimmer: Trimming barcode indexes from ${i}.${fileextension} using command:\n\t$cmd1" | tee -a $dated_log $commands_log
        fastx_clipper -h | grep 'FASTX' - 2>&1 | tee -a $dated_log $commands_log
        $cmd1  2>&1 | tee -a $dated_log
        
        #2) tagdust
        cmd2="tagdust -q -f 0.001 -s -a "$opdpath"2_"$i"_artifact.fastq -o "$opdpath$cleanfile" "$solexa_primer_adapter" "$j
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "Tagdust: Removing adapter and primer sequences from $trimfile to make $cleanfile using command:\n\t$cmd2" | tee -a $dated_log $commands_log
        $cmd2  2>&1 | tee -a $dated_log
    
    # The more typical alignOnly mode
    elif [ $trimoff == "called" ] && [ $alignonly == "called" ]
    then
        echo "Running in alignOnly mode."
        
        cmd2="tagdust -q -f 0.001 -s -a "$opdpath"2_"$i"_artifact.fastq -o "$opdpath$cleanfile" "$solexa_primer_adapter" "$j
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "Tagdust: Removing adapter and primer sequences from $trimfile to make $cleanfile using command:\n\t$cmd2" | tee -a $dated_log $commands_log
        $cmd2  2>&1 | tee -a $dated_log        
    fi
    
    


    #3) fastqc_report
    if [[ $qualityoff == "notcalled" ]]
    then
        cmd_mkdir="mkdir "$opdpath"3_"${i}"_fastqc_opd"
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "  mkdir: output directory for FastQC results\n\t$cmd_mkdir" | tee -a $dated_log $commands_log
        $cmd_mkdir  2>&1 | tee -a $dated_log
        
        cmd3="fastqc -o "$opdpath"3_"$i"_fastqc_opd --noextract "$opdpath$cleanfile
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo -e "FastQC:  Quality control from $cleanfile analyzed using command below. Output ${i}_fastqc_opd directory.\n\t$cmd3" | tee -a $dated_log $commands_log
        fastqc --version 2>&1 | tee -a $dated_log
        $cmd3  2>&1 | tee -a $dated_log
    fi
    
    
    #4) bowtie
    cmd4="bowtie -q -S --nomaqround -m 1 -p $parallel --best --chunkmbs 200 --seed 123 $bowtie1path $opdpath$cleanfile $opdpath$samfile"
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "Bowtie: Aligning $cleanfile to the genome using command below. Output $samfile\n\t$cmd4" | tee -a $dated_log $commands_log
    $cmd4  2>&1 | tee -a $dated_log $opdpath$metrics
    
    ##4) bowtie2
    #printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    #echo "  Bowtie2: Aligning $cleanfile to the genome using command below. Output $samfile" | tee -a $dated_log $commands_log
    #cmd4="bowtie2 -x "$bowtie2path" -U "$opdpath$cleanfile" -S "$opdpath$samfile" -q --phred33 --sensitive --end-to-end --seed 123 --un-gz "$opd_path$unaligned
    ##Dan's command line:  bowtie2 -x dm3 -U run416.YW-Faire-8.19_ACTTGA_L001_R1.fastq -S yw_FAIRE_rep819_hits.sam -q --phred33 --sensitive --end-to-end --seed 123
    #printf "\t" 2>&1 | tee -a $dated_log $commands_log
    #echo $cmd4 2>&1 | tee -a $dated_log $commands_log 
    #$cmd4  2>&1 | tee -a $dated_log $opdpath$metrics
#    #5) samtools compress

    cmd5="samtools view -bS -o "$opdpath$bamfile" "$opdpath$samfile
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "Samtools:  Compressing $samfile into $bamfile using command:\n\t$cmd5" | tee -a $dated_log $commands_log
    $cmd5  2>&1 | tee -a $dated_log

    #6) samtools sort
    cmd6="samtools sort -o "$opdpath$bam_sort_file" "$opdpath$bamfile 
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "Samtools:  Sorting $bamfile into $bam_sort_file using command:\n\t$cmd6" | tee -a $dated_log $commands_log
    $cmd6  2>&1 | tee -a $dated_log
    
    #7) bedtools, convert .bam -> .bed
    cmd7="bedtools bamtobed -i "$opdpath$bam_sort_file" > "$opdpath$bed_file
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "Bedtools:  Converting $bam_sort_file to $bed_file using command:\n\t$cmd7" | tee -a $dated_log $commands_log
    bedtools --version 2>&1 | tee -a $dated_log $commands_log
    bedtools bamtobed -i $opdpath$bam_sort_file > $opdpath$bed_file 2>&1 | tee -a $dated_log
    
    #8) bedToBw.sh: convert .bed -> .bw
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    inputbedfile="${i}_opd/7_${i}.bed"
    cmd8="bash bin/bedToBw.sh $inputbedfile $extension $chromlength -n -bw"
    echo -e "bedToBw:  Converting $bed_file to $bw_file using command:\n\t$cmd8" | tee -a $dated_log $commands_log
    eval $cmd8 2>&1 | tee -a $dated_log
    mv "${i}_opd/7_${i}x${extension}n.bw" "${i}_opd/8_${i}x${extension}n.bw" 2>&1 | tee -a $dated_log
    
done



#Printout a summary
printf "\n\n"  | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log
printf $(date +"%Y-%m-%d_%H:%M")"\tSUMMARY\n" | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log



#Start collecting all the data with this header: 
echo -e "Name\t#Rawreads\t#cleanReads\tPercentCleanReads\t#bowtieMappedReads\tMappedPercent\t#BowtieFailedReads\tFailedPercent\t#BowtieSuppressed\tSuppressedPercent"  | tee -a $dated_log $commands_log


#Loop through and quantify summarizing stats:
for j in ${input_files[@]}
do
    
    #Parse input file:
    path=$(dirname "${j}")
    
    myfile=$(basename $j)
    
    i="${myfile%%.*}"

    # Get the numer of rawreads in the original .fastq file:
    rawreads=($(wc $j))
    rawnum=$((rawreads / 4))
    
    
    # Read into the bowtiemetrics outputfile for bowtie stats
    opdpath="${i}_opd/"
    metrics="4_${i}_bowtiemetrics.txt"

    cleanwc=0
    cleannum=0
    cleanpercent=0
    
    # Count the nubmer of lines in the clean.fastq file if it was generated. Divide that number by 4 to get the number of reads that passed tagdust.
    if [ -e ${i}_opd/2_${i}_clean.fastq ]
    then
        cleanwc=($(wc ${i}_opd/2_${i}_clean.fastq))
        cleannum=$((cleanwc / 4))
    fi
    
    # Get the Tagdust stats.
    if [[ $rawnum > 0 ]]
    then
        cleanpercent=$((100 * cleannum / rawnum))
    fi
    
    bowtiealigned=0
    bowtiealignedpercent=0
    bowtiesuppressed="NA"
    bowtiefailed="NA"
    bowtiefailedpercent="NA"
    
    # Get the Bowtie Stats:
    if [ -e ${i}_opd/4_${i}_bowtiemetrics.txt ]
    then
        #printf "running bowtiemetrics"
        bowtiealigned=$(grep -Po "# reads with at least one reported alignment: [0-9]+" ${opdpath}${metrics} | grep -Po "[0-9]+" )
        bowtiealignedpercent=$(grep -Po "# reads with at least one reported alignment: [0-9]+ \([0-9]+.[0-9]+%\)" ${opdpath}${metrics} | grep -Po "[0-9]+.[0-9]+%" )
        bowtiefailed=$(grep -Po "# reads that failed to align: [0-9]+" ${opdpath}${metrics} | grep -Po "[0-9]+" )
        bowtiefailedpercent=$(grep -Po "# reads that failed to align: [0-9]+ \([0-9]+.[0-9]+%\)" ${opdpath}${metrics} | grep -Po "[0-9]+.[0-9]+%" )
        bowtiesuppressed=$(grep -Po "# reads with alignments suppressed due to -m: [0-9]+ \([0-9]+.[0-9]+%\)" ${opdpath}${metrics} | grep -Po "[0-9]+ " )
        bowtiesuppressedpercent=$(grep -Po "# reads with alignments suppressed due to -m: [0-9]+ \([0-9]+.[0-9]+%\)" ${opdpath}${metrics} | grep -Po "[0-9]+.[0-9]+%" )
    fi
    
    # print the output
    echo -e "${i}\t${rawnum}\t${cleannum}\t${cleanpercent}%\t${bowtiealigned}\t${bowtiealignedpercent}\t${bowtiefailed}\t${bowtiefailedpercent}\t${bowtiesuppressed}\t${bowtiesuppressedpercent}" | tee -a $dated_log $commands_log
   
done


##############################
#If cleanupOff is called, exit
##############################
if [[ $cleanoff == called ]]
then
    exit

fi
    
    
    
##############################
#CLEANUP LOOP
##############################

#Print report
printf "\n\n"  | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log
printf $(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
echo "CLEANING UP:" | tee -a $dated_log $commands_log
echo "######################################################################"  | tee -a $dated_log $commands_log
    
#for each file, compress or remove a bunch of stuff.
for i in ${list[@]}
do
    printf "\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo -e "Cleaning up ${i}" | tee -a $dated_log $commands_log
    
    #Get some file names
    opdpath="${i}_opd/"

    cleanfile="2_${i}_clean.fastq"
    artifact_file="2_${i}_artifact.fastq"
    samfile="4_${i}_output.sam"
    unaligned="4_${i}_unaligned.fastq"
    metrics="4_${i}_bowtiemetrics.txt"
    bamfile="5_${i}.bam"
    bam_sorted="6_${i}_sorted"
    bam_sort_file="6_${i}_sorted.bam"
    bed_file="7_${i}.bed"

    #gzip split fastq file if it was created
    if [[ $alignonly == "notcalled" ]]
    then
    gzip ${i}.fastq 2>&1 | tee -a $dated_log $commands_log
    fi
    
    #gzip trimfile if it was created
    if [[ $trimOff == "notcalled" ]]
    then
        trimfile="1_${i}_trim.fastq"    
        gzip $opdpath$trimfile 2>&1 | tee -a $dated_log $commands_log
    fi
    
    #gzil artifact file and bamfile
    #gzip $opdpath$cleanfile 2>&1 | tee -a $dated_log $commands_log
    gzip $opdpath$artifact_file 2>&1 | tee -a $dated_log $commands_log
    gzip $opdpath$bam_sort_file 2>&1 | tee -a $dated_log $commands_log
    
    #delete samfile, bamfile, befile, and cleanfile
    rm $opdpath/$samfile 2>&1 | tee -a $dated_log $commands_log
    rm $opdpath/$bamfile 2>&1 | tee -a $dated_log $commands_log
    rm $opdpath/$bed_file 2>&1 | tee -a $dated_log $commands_log
    rm $opdpath$cleanfile 2>&1 | tee -a $dated_log $commands_log
    
    #move
    #mv ${i}.fastq.gz ${opdpath} 2>&1 | tee -a $dated_log $commands_log

done

if [[ $alignonly == "notcalled" ]]
    then
    # rename the unmatched file
    mv unmatched.fastq ${list[0]}_unmatched.fastq 2>&1 | tee -a $dated_log $commands_log
fi