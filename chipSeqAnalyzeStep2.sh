################################################
#  ChipSeqAnalyzeStep2.sh
#  Copyright (c) 2015, Erin Osborne Nishimura
#
#PROGRAM
#       ChipSeqAnalyzeStep2.sh - To automate the peak calling and downstream formatting and analysis of multiplexed ChIP-seq data
#
#
#USAGE
#       bash ChipSeqAnalyzeStep2.sh [options] --guide <guidefile.txt>
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
#                 #name	rep	antibody    input   color
#                 EO038	REP1	ELT-2   EO042   blue
#                 EO039	REP1	IgG   EO042  red
#
#
#OPTIONS
#   --guide <guidefile.txt>     See "THE GUIDE FILE" above
#   --extension <n>             Specify the length bp to extend reads in the .wig file
#   --cleanOff                  Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
#                               It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
#   --scale <n>                 UCSC genome browser track information will include a line to turn autoscaling OFF and to set the scale value to <n>
#   -p <n>                      Runs java-genomics-toolkit in parallel. Values accepted are 1 - 8. Suggest 2 - 4. Make sure to match this number with bsub -n <n> -R "span[hosts=1]". Default=1
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
   *--guide <guidefile.txt>          
   *--extension <n>             Specify the length bp to extend reads in the .wig file. default = 150
    --cleanOff                  Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the _clean.fastq, .sam, .bam, and .bed files.
                                It retains the split.fastq.gz, _trim.fastq.gz, _sorted.bam.gz, and .bw files
    --scale <n>                 UCSC genome browser track information will include a line to turn autoscaling OFF and to set the scale value to <n>. default = 30
    -p <n>                      Runs java-genomics-toolkit in parallel. Values accepted are 1 - 8. Suggest 2 - 4. Make sure to match this number with bsub -n <n> -R \"span[hosts=1]\". Default=1
    --chromlength </path/file>  Path to the ce_chrom_length.txt file. default = /proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt
    --uploadpath </path/>       Path to an optional upload file. default = /proj/lieblab/genome_browser/ErinUCSC/
"






##################################################################################################
##################################################################################################

##########################
#1 PRE-PROCESSING: Get input and options
##########################


    #Start logfiles
        DATE=$(date +"%Y-%m-%d_%H%M")
        dated_log=${DATE}_peaks_and_tracks.log
        dated_url_log=${DATE}_url_upload.log
        dated_track_log=${DATE}_track_upload.log
        
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
        module list  2>&1 | tee -a $dated_log
       # macs2 --version
       # bash toolRunner.sh --version
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
        extension=150
        scale=30
        celength="/proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt"
        parallel=1
        uppath="/proj/lieblab/genome_browser/ErinUCSC/"
        
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
            -p) shift;
            parallel=${1:--}
            shift
            ;;
            --chromlength) shift;
            celength=${1:--}
            shift
            ;;
            --uploadpath) shift;
            uppath=${1:--}
            shift
            ;;
        esac
    done
    

        #Print out the intput files and parameters
        echo -e $(date +"%Y-%m-%d_%H:%M")"\tOPTIONS CAPTURED" | tee -a $dated_log 
        echo -e "\tGuidefile is:\t\t\t${guidefile}" | tee -a $dated_log 
        echo -e "\tBase pair extension is:\t\t${extension}" | tee -a $dated_log 
        echo -e "\tUCSC Track Scale is:\t\t${scale}" | tee -a $dated_log
        echo -e "\tChromosome Length file is:\t\t${celength}" | tee -a $dated_log
    
        #exit if required input parameters are missing
        if [ -z "$guidefile" ]
        then
            echo "ERROR: No options or input files supplied:\n" | tee -a $dated_log 
            echo "$usage"  | tee -a $dated_log 
            exit
        elif [ -z "$scale" ]
        then
            echo "ERROR: No options or input files supplied:\n" | tee -a $dated_log 
            echo "$usage"  | tee -a $dated_log 
            exit
        fi

    #Skipping commented lines, capture the number of replicates in @reparray
    #Capture the types of antibodies per replicate in @samplearray'
    #Capture the types of samples in @typearray
        filearray=($(awk 'NF &&  $1!~/^#/{print $1}' $guidefile | uniq ))
        reparray=($(awk 'NF &&  $1!~/^#/{print $2}' $guidefile ))
        uniqreparray=($(awk 'NF &&  $1!~/^#/{print $2}' $guidefile | uniq ))
        samplearray=($(awk 'NF &&  $1!~/^#/{print $3}' $guidefile ))
        inputarray=($(awk 'NF &&  $1!~/^#/{print $4}' $guidefile ))
        uniqinputarray=($(awk 'NF &&  $1!~/^#/{print $4}' $guidefile | uniq ))
        color=($(awk 'NF &&  $1!~/^#/{print $5}' $guidefile ))
        
        echo -e "\n"$(date +"%Y-%m-%d_%H:%M")"Files to process:" | tee -a $dated_log
        echo -e "\tSorted bam files for peak finding: ${filearray[*]}" | tee -a $dated_log 
        echo -e "\tReplicates: ${reparray[*]}" | tee -a $dated_log 
        echo -e "\tUnique Replicates: ${uniqreparray[*]}" | tee -a $dated_log 
        echo -e "\tSamples: ${samplearray[*]}" | tee -a $dated_log 
        echo -e "\tInput samples: ${inputarray[*]}" | tee -a $dated_log
        echo -e "\tUnique Inputs: ${uniqinputarray[*]}" | tee -a $dated_log


##########################
#2 Run MACS2
##########################
    echo -e "\n\n######################################################################"  | tee -a $dated_log
    printf "PEAK CALLING\n" | tee -a $dated_log
    echo "######################################################################"  | tee -a $dated_log
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tMACS2 initiated:\n" | tee -a $dated_log 
    
    #Set a counter
    m=0
    
    #iterate over the counter for each file in the filearray:
    while [ ${m} -lt ${#filearray[@]} ]; do
        
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tMACS2 analysis of samples:\n" | tee -a $dated_log
        #Figure out the chip_bam and input_bam paths and files:
        chip_bam=${filearray[$m]}_opd/6_${filearray[$m]}_sorted.bam
        input_bam=${inputarray[$m]}_opd/6_${inputarray[$m]}_sorted.bam
        
            #unzip if necessary
            if [ -e ${chip_bam}.gz ]; then
                echo -e "$(date +"%Y-%m-%d_%H:%M")\tgunzip: unzipping ${chip_bam}.gz" | tee -a $dated_log
                gunzip ${chip_bam}.gz  2>&1 | tee -a $dated_log
            fi
            if [ -e ${input_bam}.gz ]; then
                echo -e "$(date +"%Y-%m-%d_%H:%M")\tgunzip: unzipping ${input_bam}.gz" | tee -a $dated_log
                gunzip ${input_bam}.gz  2>&1 | tee -a $dated_log
            fi
        
        echo -e "\t\tchip_bam file is:\t${chip_bam}"
        echo -e "\t\tinput_bam file is:\t${input_bam}"
        
        #Call MACS2:
        
        quality=0.005
        cmd1="macs2 -t ${chip_bam} -c ${input_bam} -f BAM -n ${filearray[$m]}_opd/9_${filearray[$m]}_${inputarray[$m]}_${quality} -g ce --nomodel -q ${quality}"
        echo -e "${cmd1}\n" | tee -a $dated_log
        $cmd1 2>&1 | tee -a $dated_log
        
        #Convert MACS2_peak.bed files to .bam files
        cmd100="bash /proj/dllab/Erin/executables/101_bedToIndexBam/bedToIndexedBam.sh ${filearray[$m]}_opd/9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks.bed"
        echo -e "${cmd100}"  | tee -a $dated_log
        $cmd100 2>&1 | tee -a $dated_log
        
        #Move .bam files and .bam.bai files to upload area:
        cmd101="cp ${filearray[$m]}_opd/9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam ${uppath}"
        echo -e "${cmd101}"  | tee -a $dated_log
        $cmd101 2>&1 | tee -a $dated_log
        
        cmd101="cp ${filearray[$m]}_opd/9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam.bai ${uppath}"
        echo -e "${cmd101}\n"  | tee -a $dated_log
        $cmd101 2>&1 | tee -a $dated_log
        
        #Make a URL and trackfile entry for the MACS_peaks.bam file:
        echo -e "\tURL: Including 9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam in the url file $dated_url_log." | tee -a $dated_log
        echo -e "http://trackhubs.its.unc.edu/lieblab/ErinUCSC/9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam" | tee -a $dated_url_log $dated_log
        
        #Make a TRACKS file entry for MACS_peaks.bam:
        
        rep=${reparray[$m]}
        sample=${samplearray[$m]}
        
        echo -e "\tTracks: Including 9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam in the track file ${dated_track_log}." | tee -a $dated_log
        echo -e "track name=${sample}_${rep}_macs2 description=9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted bamcolormode=${color[$m]} bigDataUrl=http://trackhubs.its.unc.edu/lieblab/ErinUCSC/9_${filearray[$m]}_${inputarray[$m]}_${quality}_peaks_sorted.bam type=bam" | tee -a $dated_track_log $dated_log
        
        echo -e "\n\n" | tee -a $dated_log
        #Iterate counter
        let m=m+1
        
    done



##########################
#3 Modify the .bw output file to match sample types. Generate a track URL file.
##########################

    echo -e "\n\n######################################################################"  | tee -a $dated_log
    printf "GENERATING TRACK FILES\n" | tee -a $dated_log
    echo "######################################################################"  | tee -a $dated_log
    
    
    
    
    n=0
    #Loop through the chip samples, move .bw files, and generate url and track files:
    while [ ${n} -lt ${#filearray[@]} ]; do
        
        chip=${filearray[$n]}
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tGenerating track files for: ${chip}\n" | tee -a $dated_log
        
    
        ####### copy files to upload area ###########
        echo -e "\tCopying 8_${chip}x${extension}n.bw to ftp site." | tee -a $dated_log
        cp ${chip}_opd/8_${chip}x${extension}n.bw ${uppath} 2>&1 | tee -a $dated_log
        
        ####### compile a list of url's ###########
        echo -e "\tURL: Including 8_${chip}x${extension}n.bw in the url file $dated_url_log." | tee -a $dated_log
        echo -e "http://trackhubs.its.unc.edu/lieblab/ErinUCSC/8_${chip}x${extension}n.bw" | tee -a $dated_url_log $dated_log
        
        
        ####### compile a list of track names ###########
        #Get replicate #
        rep=${reparray[$n]}
        
        #Get Sample Name
        sample=${samplearray[$n]}
        
        #Get color
        if [[ ${color[$n]} == "blue" ]]; then
            color="0,0,255"
        elif [[ ${color[$n]} == "red" ]]; then
            color="255,0,0"
        elif [[ ${color[$n]} == "green" ]]; then
            color="0,255,0"
        fi
        
        echo -e "\tTracks: Including 8_${chip}x${extension}n.bw in the track file ${dated_track_log}." | tee -a $dated_log
        echo -e "track name=${sample}_${rep} desc=8_${chip}x${extension}n maxHeightPixels=100:100:11 autoScale=off viewLimits=0:${scale} color=${color} visibility=2 bigDataUrl=http://trackhubs.its.unc.edu/lieblab/ErinUCSC/8_${chip}x${extension}n.bw type=bigWig" | tee -a $dated_track_log $dated_log
        
        ((n+=1))
    done
    
    
    #Loop through the input samples, move .bw files, and generate url and track files:
    u=0
    #Loop through the input samples:
    while [ ${u} -lt ${#uniqinputarray[@]} ]; do
        
        chip=${uniqinputarray[$u]}
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tGenerating track files for: ${chip}\n" | tee -a $dated_log
        
    
        ####### copy files to upload area ###########
        echo -e "\tCopying 8_${chip}x${extension}n.bw to ftp site." | tee -a $dated_log
        cp ${chip}_opd/8_${chip}x${extension}n.bw /proj/lieblab/genome_browser/ErinUCSC/ 2>&1 | tee -a $dated_log
        
        ####### compile a list of url's ###########
        echo -e "\tURL: Including 8_${chip}x${extension}n.bw in the url file $dated_url_log." | tee -a $dated_log
        echo -e "http://trackhubs.its.unc.edu/lieblab/ErinUCSC/8_${chip}x${extension}n.bw" | tee -a $dated_url_log $dated_log
        
        
        ####### compile a list of track names ###########
        #Get replicate #
        rep=
        
        num=0
        echo -e "starting rep loop\n"
        
        for j in ${inputarray[@]}; do
            if [[ ${j} == ${uniqinputarray[$u]} ]]; then 
                rep=${num}
                break
            fi

            num=$((num + 1))
        done
        
        
        echo -e "\tTracks: Including 8_${chip}x${extension}0n.bw in the track file ${dated_track_log}." | tee -a $dated_log
        echo -e "track name=input_${reparray[$num]} desc=8_${chip}x${extension}n maxHeightPixels=100:100:11 autoScale=off viewLimits=0:${scale} color=0,0,0 visibility=2 bigDataUrl=http://trackhubs.its.unc.edu/lieblab/ErinUCSC/8_${chip}x${extension}n.bw type=bigWig" | tee -a $dated_track_log $dated_log
        
        ((u+=1))
    done


##########################
#4 Generate subtraction files. This loop calls java-genomics-toolkit.
##########################

    x=0
    
    
    while [ ${x} -lt ${#filearray[@]} ]; do
        
        chip=${filearray[$x]}
        input=${inputarray[$x]}
        igG=
        
        echo -e "filearray value is ${filearray[$x]}"
        
        y=0
        while [ ${y} -lt ${#reparray[@]} ]; do
            #echo -e "\ty value is ${y}"
            if [ ${reparray[$y]} == ${reparray[$x]} ]; then
                #echo -e "\t REPS MATCH"
                #echo -e "\t\t${samplearray[$y]}"
                if [ "${samplearray[$y]}" == "IgG" ]; then                    
                    IgG=${filearray[$y]}
                    
                fi
                
            fi
            
            ((y+=1))
        done
        
        echo -e "outside of the loop. IgG file for rep ${reparray[$x]} is ${IgG}"
            
        chip_bw=${chip}_opd/8_${chip}x${extension}n.bw
        igg_bw=${IgG}_opd/8_${IgG}x${extension}n.bw
        #
        ##JavaGenomicsToolkit:
        
        echo -e "IgG is ${IgG} and Chip is ${chip}"
        if [ ${IgG} != ${chip} ]; then
            echo -e "\n\n"$(date +"%Y-%m-%d_%H:%M")"\tjava-genomics-toolkit: Generating Subtraction bigWig files between Chip and IgG samples:" | tee -a $dated_log
            
            cmd5="bash /proj/dllab/Erin/executables/java-genomics-toolkit-2/java-genomics-toolkit/toolRunner.sh wigmath.Subtract -m ${chip_bw} -s ${igg_bw} -o ${chip}_opd/${chip}_over_${IgG}.wig -f -p ${parallel}"
            echo -e "\n${cmd5}" | tee -a $dated_log
            #$cmd5 2>&1 | tee -a $dated_log
            
            ##Convert to bw:
            cmd6="wigToBigWig ${chip}_opd/${chip}_over_${IgG}.wig ${celength} ${chip}_opd/${chip}_over_${IgG}.bw"
            echo -e "\n${cmd6}" | tee -a $dated_log
            $cmd6 2>&1 | tee -a $dated_log
            
            echo -e "\tCopying ${chip}_opd/${chip}_over_${IgG}.bw to ftp site." | tee -a $dated_log
            cp ${chip}_opd/${chip}_over_${IgG}.bw /proj/lieblab/genome_browser/ErinUCSC/ 2>&1 | tee -a $dated_log
            
            #add to url file
            echo -e "\tIncluding ${chip}_over_${IgG}.bw in the url file $dated_url_log." | tee -a $dated_log
            echo -e "http://trackhubs.its.unc.edu/lieblab/ErinUCSC/${chip}_over_${IgG}.bw" | tee -a $dated_url_log    
            
            #add to track file
            rep=${reparray[$x]}
            sample=${samplearray[$x]}
            echo -e "\tIncluding ${chip}_over_${IgG}.bw in the track file ${dated_track_log}" | tee -a $dated_log
            echo -e "track name=${sample}_${rep}_minus_${IgG} description=${chip}_over_${IgG} maxHeightPixels=100:100:11 autoScale=off viewLimits=0:$((scale / 5)) color=0,0,0 visibility=2 bigDataUrl=http://trackhubs.its.unc.edu/lieblab/ErinUCSC/${chip}_over_${IgG}.bw type=bigWig" | tee -a $dated_track_log
            #
        fi
        
        ((x+=1))
    done
    

    
    
#PREVIOUS LOOPS/SCRIPTS:
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