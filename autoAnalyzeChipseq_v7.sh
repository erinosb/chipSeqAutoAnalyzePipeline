#! /bin/sh/

################################################
#  autoAnalyzeChipseq_v5.sh
#
#PROGRAM
#   autoAnalyzeChipseq_v5.sh - To automate the analysis of multiplexed ChIP-seq data
#
#
#USAGE
#
#   bash autoAnalyzeChipseq_v5.sh [options] <inputFile.txt> <barcodeIndexFile.txt>
#   OR
#   bash autoAnalyzeChipseq_v5.sh --splitOnly [options] <inputFile.txt> <barcodeIndexFile.txt>
#   OR
#   bash autoAnalyzeChipseq_v5.sh --splitOFF [options] <input1.fastq> input2.fastq
#
#
#MODES
#   There are three modes for use... --split mode (default), --splitOnly and --splitOff mode
#   --split mode (default)      This mode uses as input a single homebrew multiplexed .txt file and
#                               a single barcode index file. The pipeline then...
#                                   0) splits the multiplexed file into single sample files by index
    #                               1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
    #                               2) perform Tagdust using a solexa library of adapter and primer sequences.
    #                               3) perform Fastqc to make a report of quality
    #                               4) bowtie alignment
    #                               5) compress .sam --> .bam using samtools view
    #                               6) Sort .bam into _sorted.bam using samtools sort
    #                               7) Convert _sorted.bam into .bed using bedtools
    #                               8) Convert .bed into .wig using zinba.
    #                               9) gzip .wig to .wig.gz
    #
#   --splitOnly                 This mode uses as input a single homebrew multiplexed .txt file and
#                               a single barcode index file. The pipeline then...
#                                   0) splits the multiplexed file into single sample files by index
#
#   --splitOff                  This mode uses as input a list of .txt or .fastq files and analyzes them...
#                               a single barcode index file. The pipeline then...
    #                               1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
    #                               2) perform Tagdust using a solexa library of adapter and primer sequences.
    #                               3) perform Fastqc to make a report of quality
    #                               4) bowtie alignment
    #                               5) compress .sam --> .bam using samtools view
    #                               6) Sort .bam into _sorted.bam using samtools sort
    #                               7) Convert _sorted.bam into .bed using bedtools
    #                               8) Convert .bed into .wig using zinba.
    #                               9) gzip .wig to .wig.gz
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
#OPTIONS
#
#   --splitOnly                 Runs in spligOnly mode. Suppresses analysis.
#   --splitOff                  Runs in splitOff mode. Suppresses splitting sequencefiles
#   --qualityOff                Suppresses quality score reports
#   --extension <n>             Specify the length bp to extend reads in the .wig file
#   --trimOff                   Suppresses trimming six basepairs of multiplexing indices
#   
#AUTHOR
#   Erin Osborne Nishimura
#
#DATE
#   April 4, 2014
#
#BUGS
#   -- check whether certain modules have been loaded;
#   -- Auto load all required modules
#
#
#TO FIX
#
#POTENTIAL FUTURE EXPANSION
#   Accept batch multiplex sequencefiles
#   Switch to accept both inputFile.txt and inputFile.fastq and to croak and die if neither suffixes are present
#   Change the optioning such that there are three different modes... --split, --analyze, and --cleanup
#
#################################################

#####################   SET VARIABLES   ######################
solexa_primer_adapter="/proj/dllab/Erin/sequences/solexa-library-seqs.fasta"    #tagdust needs a .fasta file that contains a list of all the solexa primer and adapter sequences. Set this
                                                                                  #variable to a path pointing to that file.
bowtie2path="/proj/dllab/Erin/ce10/from_ucsc/seq/genome_bt2/ce10"               #tophat needs to know where the bowtie2 index files are located. Set this varaible to the path and root
                                                                                  #name of those index files.
                                                                                  #Also, the genome sequence (a .fa file) also needs to be in that same directory.
extension=100                                                                         #zinba needs to know how long your reads are so that it can make a .wig file with the proper extension lengths
twobit=/proj/dllab/Erin/ce10/from_ucsc/seq/ce10.2bit                            #zinba needs to know where the bowtie twobit files are located. I'm not sure whether zinba works with updated bowtie2
                                                                                  #files or whether you need the old bowtie twobit files.
                                                                                  
##################################################################################################



usage="
PROGRAM
   autoAnalyzeChipseq_v5.sh - To automate the analysis of multiplexed ChIP-seq data


USAGE

   bash autoAnalyzeChipseq_v5.sh [options] <inputFile.txt> <barcodeIndexFile.txt>
   OR
   bash autoAnalyzeChipseq_v5.sh --splitOnly [options] <inputFile.txt> <barcodeIndexFile.txt>
   OR
   bash autoAnalyzeChipseq_v5.sh --splitOFF [options] <input1.fastq> input2.fastq


MODES
   There are three modes for use... --split mode (default), --splitOnly and --splitOff mode
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
                                   8) Convert .bed into .wig using zinba.
                                   9) gzip .wig to .wig.gz
    
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
                                   8) Convert .bed into .wig using zinba.
                                   9) gzip .wig to .wig.gz


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

   --splitNAlign                Runs this script in a split and align mode.
   --splitOnly                  Runs this script in splitOnly mode. Suppresses alignment.
   --alignOnly                  Runs in splitOff mode. Suppresses splitting sequencefiles
   --qualityOff                 Suppresses quality score reports
   --extension <n>              Specify the length bp to extend reads in the .wig file
   --trimOff                    Suppresses trimming six basepairs of multiplexing indices"


#################
#PRE-PROCESSING: Load modules, Check errors, get filenames, set options
#################


#Load modules
#module add samtools
#module add bedtools
#module add bowtie2
#module add fastqc



#Start logfiles
DATE=$(date +"%Y-%m-%d_%H%M")
dated_log=${DATE}.log
commands_log=${DATE}_commands.log
echo $DATE | tee -a $dated_log $commands_log

echo "INITIATED autoAnalyzeChipseq_v2.sh using command: $0 $*" | tee -a $dated_log $commands_log
printf "\n" | tee -a $dated_log $commands_log

if [ -z "$1" ]
  then
    echo "ERROR: No inputFile supplied:" | tee -a $dated_log $commands_log
    echo "$usage"  | tee -a $dated_log $commands_log
    exit
fi

if [ -z "$2" ]
  then
    echo "ERROR: No barcodeIndexFile supplied:" | tee -a $dated_log $commands_log
    echo "$usage"  | tee -a $dated_log $commands_log
    exit
fi

#Set options

alignonly="notcalled"
qualityoff="notcalled"
splitonly="notcalled"
trimoff="notcalled"
list=$@

#Get options
for n in $@
do
    case $n in
        --splitOff) shift;
        splitoff="called"
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
    esac
done



#Un-multiplex using barcode splitter
if [[ $splitoff == "notcalled" ]]
then
    multi=$1
    index=$2

    printf "\nSplitting seqfile:\n  "  | tee -a $dated_log $commands_log
    
    echo $multi | tee -a $dated_log $commands_log
    printf "into multiple files based on barcodes in:\n  " | tee -a $dated_log $commands_log
    echo $index | tee -a $dated_log $commands_log
    printf "using command:\n\n  " | tee -a $dated_log $commands_log
    
    
    
    ######################
    #Split the multiplexed file into multiple files based on barcoded indexes
    ######################
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\t"
    echo "cat $multi | fastx_barcode_splitter.pl --bcfile $index --prefix "" --suffix ".fastq" --bol" | tee -a $dated_log $commands_log
    cat $multi | fastx_barcode_splitter.pl --bcfile $index --prefix "" --suffix ".fastq" --bol | tee -a $dated_log $commands_log

    #####################
    #Stop the program here if --splitOnly is called
    #####################
    if [[ $splitonly == "called" ]]
    then
        exit
    fi
    
    
    
    ######################
    #Loop through each of the new split sequence files and perform tasks on individuals:
    ######################

    list=$(awk '{print $1;}' $index | grep -v '#')
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\t"
    printf "\nWill Chip-seq analyze the following .fastq/.txt files:\n" | tee -a $dated_log $commands_log
    for i in $list
    do
        echo ${i}\.fastq | tee -a $dated_log $commands_log
    done
    printf "\n" | tee -a $dated_log $commands_log
fi


##############################
#If splitOnly is called, exit
##############################

if [[ $splitonly == "called" ]]
then
    exit
fi

#############################
#If splitOff is called, use the list of remaining arguments as input files...
#############################

if [[ $splitoff == "called" ]]
then
    list=$@
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t\t"
    printf "\nWill Chip-seq analyze the following .fastq/.txt files:\n" | tee -a $dated_log $commands_log
    for i in $list
        do
            echo ${i} | tee -a $dated_log  $commands_log
        done
fi



#############################
#Quality Processing
#############################

#For each sample,
    #1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
    #2) perform Tagdust using a solexa library of adapter and primer sequences.
    #3) perform Fastqc to make a report of quality
    #4) bowtie alignment
    #5) compress .sam --> .bam using samtools view
    #6) Sort .bam into _sorted.bam using samtools sort
    #7) Convert _sorted.bam into .bedgraph (.wig file) using zinba
    

for j in $list
do
    #Remove the suffix
    i="${j%.*}"
    
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
    wig_file="8_${i}.wig"
    r_code="8_${i}_code.R"
    processlog=${i}_process.log
    btlog=${i}_bt.log
    

    

    printf "\n\n"  | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log
    printf $(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "PROCESSING $i:" | tee -a $dated_log $commands_log
    echo "######################################################################"  | tee -a $dated_log $commands_log
    
    #Make Output Directory
    printf "mkdir $opd" 2>&1 | tee -a $dated_log 
    mkdir $opd 2>&1 | tee -a $dated_log

    if [[ $trimoff == "notcalled" ]]
    then
        #1) fastx_trim:
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo "  fastx_trimmer: Trimming barcode indexes from ${i}.fastq using command:" | tee -a $dated_log $commands_log
        cmd1="fastx_trimmer -f 9 -Q 33 -i "$i".fastq -o "$opdpath$trimfile
        printf "\t" 2>&1 | tee -a $dated_log $commands_log
        echo $cmd1 2>&1 | tee -a $dated_log $commands_log
        fastx_clipper -h | grep 'FASTX' - 2>&1 | tee -a $dated_log 
        $cmd1  2>&1 | tee -a $dated_log
    fi

    #2) tagdust
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "  Tagdust: Removing adapter and primer sequences from $trimfile to make $cleanfile using command:" | tee -a $dated_log $commands_log
    cmd2="tagdust -q -f 0.001 -s -a "$opdpath"2_"$i"_artifact.txt -o "$opdpath$cleanfile" "$solexa_primer_adapter" "$opdpath$trimfile
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    echo $cmd2 2>&1 | tee -a $dated_log $commands_log
    $cmd2  2>&1 | tee -a $dated_log

    #3) fastqc_report
    if [[ $qualityoff == "notcalled" ]]
    then
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        
        cmd_mkdir="mkdir "$opdpath"3_"${i}"_fastqc_opd"
        echo $cmd_mkdir 2>&1 | tee -a $dated_log $commands_log
        $cmd_mkdir  2>&1 | tee -a $dated_log
        
        printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
        echo "  Fastqc:  Quality control from $cleanfile analyzed using command below. Output $i _fastqc_opd directory." | tee -a $dated_log $commands_log
        fastqc --version 2>&1 | tee -a $dated_log
        cmd3="fastqc -o "$opdpath"3_"$i"_fastqc_opd --noextract "$opdpath$cleanfile
        printf "\t" 2>&1 | tee -a $dated_log $commands_log
        echo $cmd3 2>&1 | tee -a $dated_log $commands_log
        $cmd3  2>&1 | tee -a $dated_log
    fi
    
    
    #4) bowtie
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "  Bowtie2: Aligning $cleanfile to the genome using command below. Output $samfile" | tee -a $dated_log $commands_log
    cmd4="bowtie2 -x "$bowtie2path" -U "$opdpath$cleanfile" -S "$opdpath$samfile" -q --phred33 --sensitive --end-to-end --seed 123 --un-gz "$opd_path$unaligned
    #Dan's command line:  bowtie2 -x dm3 -U run416.YW-Faire-8.19_ACTTGA_L001_R1.fastq -S yw_FAIRE_rep819_hits.sam -q --phred33 --sensitive --end-to-end --seed 123
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    echo $cmd4 2>&1 | tee -a $dated_log $commands_log 
    $cmd4  2>&1 | tee -a $dated_log $opdpath$metrics
    
    #5) samtools compress
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "  Samtools:  Compressing $samfile into $bamfile using command:" | tee -a $dated_log $commands_log
    cmd5="samtools view -bS -o "$opdpath$bamfile" "$opdpath$samfile
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    echo $cmd5 2>&1 | tee -a $dated_log $commands_log
    $cmd5  2>&1 | tee -a $dated_log

    #6) samtools sort
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "  Samtools:  Sorting $bamfile into $bam_sort_file using command:" | tee -a $dated_log $commands_log
    cmd6="samtools sort "$opdpath$bamfile" "$opdpath$bam_sorted
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    echo $cmd6 2>&1 | tee -a $dated_log $commands_log
    $cmd6  2>&1 | tee -a $dated_log
    
    #7) bedtools, convert .bam -> .bed
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "   Bedtools:  Converting $bam_sort_file to $bed_file using command:" | tee -a $dated_log $commands_log
    cmd7="bedtools bamtobed -i "$opdpath$bam_sort_file" > "$opdpath$bed_file

    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    echo $cmd7 2>&1 | tee -a $dated_log $commands_log
    bedtools --version 2>&1 | tee -a $dated_log 
    bedtools bamtobed -i $opdpath$bam_sort_file > $opdpath$bed_file 2>&1 | tee -a $dated_log
    
    
    
    #8) zinba .bed -> .wig
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    echo "  R: Writing a .R file to convert $bed_file to $wig_file using Zinba(basealigncounts).\n" | tee -a $dated_log $commands_log
    
    zinba_code="
library(zinba)
basealigncount(
    inputfile=\"$opdpath$bed_file\",
    outputfile=\"$opdpath$wig_file\",
    extension="$extension",
    filetype=\"bed\",
    twoBitFile=\""$twobit"\"
)"
    echo "$zinba_code" > $opdpath$r_code 2>&1 | tee -a $dated_log $commands_log
    cmd8="/proj/.test/roach/FAIRE/bin/R --vanilla < "$opdpath$r_code
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    printf "%s" $cmd8 2>&1 | tee -a $dated_log $commands_log
    /proj/.test/roach/FAIRE/bin/R --vanilla < $opdpath$r_code 2>&1 | tee -a $dated_log
    printf "\n" | tee -a $dated_log
    
    #9) gzip .wig -> .wig.gz
    printf "\n\n"$(date +"%Y-%m-%d_%H:%M")"\t" | tee -a $dated_log $commands_log
    printf "  gzip: Compressing $wig_file to $wig_file.gz using command:\n" | tee -a $dated_log $commands_log
    printf "\t" 2>&1 | tee -a $dated_log $commands_log
    printf "gzip "$opdpath$wig_file 2>&1 | tee -a $dated.log $commands_log
    gzip $opdpath$wig_file 2>&1 | tee -a $dated.log $commands_log

    
done

