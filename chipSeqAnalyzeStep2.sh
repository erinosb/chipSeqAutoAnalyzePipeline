################################################
#  ChipSeqAnalyzeStep2.sh
#  Copyright (c) 2015, Erin Osborne Nishimura
#
#PROGRAM
#       ChipSeqAnalyzeStep2.sh - To automate the peak calling and downstream formatting and analysis of multiplexed ChIP-seq data
#
#
#USAGE
#       bash ChipSeqAnalyzeStep2.sh [options] --guide <guidefile.txt> -p <n>
#
#
#MODES
#
#
#
#ARGUMENTS
#       <guidefile.txt>
#
#
#
#THE GUIDE FILE
#       The Guide file is a tab deliminated text file that contains information regarding each sample within a multiplexed set of chip-seq experiments:
#                   The first line is commented with # sign and contains the titles of each column: name, replicate, antibody and type
#                   Name is the root name of each sample as specified in the barcode index file from ChipSeqAnalyzeStep1.sh
#                   Replicate is the replicate name and it must be in the form: REP#
#                   Antibody is a string describing the antibody used in the Chip-seq experiment or the word 'input' if the sample is an input sample
#                   Type is 'c' for Chip-seq samples and 'i' for input samples
#
#       Example guidefile.txt:
#
#                 #name	rep	antibody    type    color	
#                 EO038	REP1	ELT-2   c   blue
#                 EO039	REP1	IgG c   red	
#                 EO042	REP1	input   i   black
#
#
#OPTIONS
#   --guide <guidefile.txt>          
#   --extension <n>             Specify the length bp to extend reads in the .wig file
#   --cleanOff                  Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
#                               It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
#   --scale <n>                 UCSC genome browser track information will include a line to turn autoscaling OFF and to set the scale value to <n>
#   
#   
#AUTHOR
#   Erin Osborne Nishimura
#
#DATE
#   April 20, 2015
#
#BUGS/FUTURE EXPANSION
#   -- check whether certain modules have been loaded
#   -- Auto load all required modules
#   -- toggle between bowtie or bowtie2
#
#

usage="  ChipSeqAnalyzeStep2.sh
  Copyright (c) 2015, Erin Osborne Nishimura

PROGRAM
       ChipSeqAnalyzeStep2.sh - To automate the peak calling and downstream formatting and analysis of multiplexed ChIP-seq data


USAGE
       bash ChipSeqAnalyzeStep2.sh [options] --guide <guidefile.txt> -p <n>

REQUIRED ARGUMENTS
       <guidefile.txt>
       
THE GUIDE FILE
       The Guide file is a tab deliminated text file that contains information regarding each sample within a multiplexed set of chip-seq experiments:
                   The first line is commented with # sign and contains the titles of each column: name, replicate, antibody and type
                   Name is the root name of each sample as specified in the barcode index file from ChipSeqAnalyzeStep1.sh
                   Replicate is the replicate name and it must be in the form: REP#
                   Antibody is a string describing the antibody used in the Chip-seq experiment or the word 'input' if the sample is an input sample
                   Type is 'c' for Chip-seq samples and 'i' for input samples
                   Color is the color you want the track to appear in UCSC genome browser

       Example guidefile.txt:
                 #name	rep	antibody    type    color	
                 EO038	REP1	ELT-2   c   blue
                 EO039	REP1	IgG c   red	
                 EO042	REP1	input   i   black


OPTIONS
   --guide <guidefile.txt>          
   --extension <n>             Specify the length bp to extend reads in the .wig file
   --cleanOff                  Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
                               It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
   --scale <n>                 UCSC genome browser track information will include a line to turn autoscaling OFF and to set the scale value to <n>"






##################################################################################################
##################################################################################################

##########################
#1 PRE-PROCESSING: Get input and options
##########################


    #Start logfiles
        DATE=$(date +"%Y-%m-%d_%H%M")
        dated_log=${DATE}_peaks_and_tracks.log
        
    #Report that the program is running. Report $0.
        echo "######################################################################"  | tee -a $dated_log
        printf "INITIATED\n" | tee -a $dated_log
        echo "######################################################################"  | tee -a $dated_log
        printf $(date +"%Y-%m-%d_%H:%M")"\tInitiated with the following command:\n\t\t$0 $*\n" | tee -a $dated_log 
        
    #Make sure the proper modules are uploaded:
        echo -e "\n\n######################################################################"  | tee -a $dated_log
        printf "PRE-PROCESSING\n" | tee -a $dated_log
        echo "######################################################################"  | tee -a $dated_log
        printf $(date +"%Y-%m-%d_%H:%M")"\tRequired modules: samtools, bedtools, macs2" | tee -a $dated_log  
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tThis pipeline will run with the following modules:\n" | tee -a $dated_log 
        source /nas02/apps/Modules/default/init/bash
        module list  | tee -a $dated_log
       # macs2 --version
        printf "\n\n" | tee -a $dated_log

    #Warn or die if specified options are not present:
    #Check that there is some input:
        if [ -z "$1" ]
        then
            echo "ERROR: No options or input files supplied:" | tee -a $dated_log 
            echo "$usage"  | tee -a $dated_log 
            exit
        fi
    
    #Capture options
        guidefile=
        extension=
        scale=
        
        for n in $@
        do
            case $n in
            --guide) shift;
            guidefile="$1"
            shift
            ;;
            --extension) shift;
            extension=${1:--}
            shift
            ;;
            --scale) shift;
            scale=${1:--}
            shift
            ;;
        esac
    done
        
        echo -e $(date +"%Y-%m-%d_%H:%M")"OPTIONS CAPTURED" | tee -a $dated_log 
        echo -e "\tGuidefile is:\t\t\t${guidefile}" | tee -a $dated_log 
        echo -e "\tBase pair extension is:\t\t${extension}" | tee -a $dated_log 
        echo -e "\tUCSC Track Scale is:\t\t${scale}" | tee -a $dated_log 


    #Skipping commented lines, capture the number of replicates in @reparray
    #Capture the types of antibodies per replicate in @samplearray'
    #Capture the types of samples in @typearray
        filearray=($(awk 'NF &&  $1!~/^#/{print $1}' $guidefile | uniq ))
        reparray=($(awk 'NF &&  $1!~/^#/{print $2}' $guidefile ))
        uniqreparray=($(awk 'NF &&  $1!~/^#/{print $2}' $guidefile | uniq ))
        samplearray=($(awk 'NF &&  $1!~/^#@/{print $3}' $guidefile ))
        inputarray=($(awk 'NF &&  $1!~/^#/{print $4}' $guidefile ))
        color=($(awk 'NF &&  $1!~/^#/{print $5}' $guidefile ))
        
        echo -e "will process the following files: ${filearray[*]}"
        echo -e "reparray is: ${reparray[*]}"
        echo -e "uniqreparray is: ${uniqreparray[*]}"
        echo -e "samplearray is: ${samplearray[*]}"
        echo -e "inputarray is: ${inputarray[*]}"


##########################
#2 Run MACS2
##########################
    echo -e "\n\n######################################################################"  | tee -a $dated_log
    printf "PEAK CALLING\n" | tee -a $dated_log
    echo "######################################################################"  | tee -a $dated_log
    

    #Set a counter
    m=0
    
    #iterate over the counter for each file in the filearray:
    while [ ${m} -lt ${#filearray[@]} ]; do
        
        #Figure out the chip_bam and input_bam paths and files:
        chip_bam=${filearray[$m]}_opd/6_${filearray[$m]}_sorted.bam
        input_bam=${inputarray[$m]}_opd/6_${inputarray[$m]}_sorted.bam
        
            #unzip if necessary
            if [ -z ${chip_bam}.gz ]; then
                gunzip ${chip_bam}.gz
            fi
            if [ -z ${input_bam}.gz ]; then
                gunzip ${input_bam}.gz
            fi
        
        #Call MACS2:
        
        echo -e "chip_bam is ${chip_bam}"
        echo -e "inputbam is ${input_bam}"
        
        q=0.01
        cmd1="macs2 -t ${chip_bam} -c ${input_bam} -f BAM -n ${filearray[$m]}_${inputarray[$m]}_${q} -g ce --nomodel --shiftsize=125 -q ${q}"
        cmd1
        let m=m+1 
    done
    
    
    

##############################################################################
#Author: Adam Robinson
#Date: 2015-04-01
#Purpose: run sorted .bam files against input control .bam files through macs2
##############################################################################
####step 1: create array containing all files listed in elements_list.txt###
#elements=$1
#ChIP=($(cut -f1 $elements))
#Input=($(cut -f2 $elements))
#index=1
#count=$(( ${#ChIP[@]}))
#
##step 2-3 performed through while loop#
####step 2: find .bam file; find control .bam file, the last file; the control .bam file should be the last file in the options###
####step 3: call macs2 using sorted.bam file and control .bam file###
#
#while [ "$index" -lt "$count" ]
#do
#        bamChIP=$(find -name "6_${ChIP[$index]}_sorted.bam")
#        bamInput=$(find -name "6_${Input[$index]}_sorted.bam")
#        macs2 -t $bamChIP -c $bamInput -f BAM -n ${ChIP[$index]}_${Input[$index]}_q1 -g ce --nomodel --shiftsize=125
#        index=$(( $index + 1 ))
#done



##########################
#3 Modify the .bw output file to match sample types. Generate a track URL file.
##########################



