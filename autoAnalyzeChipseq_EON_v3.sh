#! /bin/bash/

################################################
#unmultiplex_EON_v1.sh
#
#PROGRAM
#   autoAnalyzeChipseq_v2.sh - To automate the analysis of multiplexed ChIP-seq data
#
#USAGE
#   bash autoAnalyzeChipseq_v2.sh <inputFile.txt> <barcodeIndexFile.txt> [options]
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
#
#OPTIONS
#
#   --splitOff                  Don't split sequencefiles
#   --qualityOff                Don't output quality score reports
#   
#AUTHOR
#   Erin Osborne Nishimura
#
#DATE
#   January 5, 2013
#
#BUGS
#
#
#TO FIX
#
#
#POTENTIAL FUTURE EXPANSION
#   Accept batch multiplex sequencefiles
#   Switch to accept both inputFile.txt and inputFile.fastq and to croak and die if neither suffixes are present
#   Include help menu
#   Include quality control
#
#################################################

usage="
    USAGE
        bash autoAnalyzeChipseq_v2.sh <inputFile.txt> <barcodeIndexFile.txt> [options]
        
    ARGUMENTS
        <inputFile.txt>           This is the file from the sequencing facility. It is a fastq file containing Illumina sequencing reads generated from multiplexed samples
        <barcodeIndexFile.txt>    This is a file listing the barcodes that are at the 5' end of each sequence string.  It should be in the format:

                   #barcode file for the multiplexed sequences generated 11/11/12
                   EO33	AGATGGT
                   EO34	CACGTCG
                   EO35	GATCTTG
                   EO36	TCAGGAC
                   EO37	ACAGTTG

    OPTIONS
           --qualityOff                Don't output quality score reports"



#################
#PRE-PROCESSING: Check errors, get filenames, set options
#################

echo "INITIATED autoAnalyzeChipseq_v2.sh using command: $0"
printf "\n"

if [ -z "$1" ]
  then
    echo "ERROR: No inputFile supplied:"
    echo "$usage"
    exit
fi

if [ -z "$2" ]
  then
    echo "ERROR: No barcodeIndexFile supplied:"
    echo "$usage"
    exit
fi

multi=$1
index=$2
multi_root=${MULTI%.txt}
#MULTI_root=${MULTI%.fastq} <-- This doesn't work. Change the program here to accept .fastq or .txt


printf "\nSplitting seqfile:\n  "
echo $multi
printf "into multiple files based on barcodes in:\n  "
echo $index
printf "using command:\n\n  "



######################
#Split the multiplexed file into multiple files based on barcoded indexes
######################

#echo "cat $multi | fastx_barcode_splitter.pl --bcfile $index  --prefix "" --suffix ".fastq" --bol"
#cat $multi | fastx_barcode_splitter.pl --bcfile $index --prefix "" --suffix ".fastq" --bol 1>> split.log 2>&1
#%%%%%%%%%%%%
#TURN THE CAT FUNCTION BACK ON!!!
#%%%%%%%%%%%%%

######################
#Quality Control
######################


#Loop through each of the new split sequence files and perform tasks on individuals:
list=$(awk '{print $1;}' $index | grep -v '#')

printf "\nWill clean up the following .fastq files:\n"
for i in $list
do
    echo ${i}\.fastq
done
printf "\n"


#For each sample,
    #1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
    #2) perform Tagdust using a solexa library of adapter and primer sequences.
    #3
    #4) bowtie alignment
    #5) compress .sam --> .bam using samtools view
    #6) Sort .bam into _sorted.bam using samtools sort
    #7) Convert _sorted.bam into .bedgraph (.wig file) using bedtools genomecov
    

for i in $list
do
    trimfile=${i}_trim.fastq
    cleanfile=${i}_clean.fastq
    samfile=${i}_output.sam
    bamfile=${i}.bam
    bam_sorted=${i}_sorted
    bam_sort_file=${i}_sorted.bam
    bed_file=${i}.bed
    processlog=${i}_process.log
    btlog=${i}_bt.log
    
    echo "PROCESSING $i:"
    
    #fastx_trim:
    echo "  fastx_trimmer: Trimming barcode indexes from ${i}.fastq using command:"
    echo "    fastx_trimmer -f 9 -Q 33 -i ${i}.fastq -o $trimfile"
    fastx_trimmer -f 9 -Q 33 -i ${i}.fastq -o $trimfile
    

    #tagdust
    echo "  Tagdust: Removing adapter and primer sequences from $trimfile using command:"
    echo "    tagdust -q -f 0.001 -a ${i}_artifact.txt -o $cleanfile /proj/dllab/solexa-library-seqs.fasta $trimfile"
    tagdust -q -f 0.001 -a ${i}_artifact.txt -o $cleanfile /proj/dllab/solexa-library-seqs.fasta $trimfile 1>> process.log 2>&1
    
    #fastqc_report
    mkdir ${i}_fastqc_opd
    echo "  Fastqc:  Quality control from ${i} analyzed using command:"
    echo "    fastqc -o ${i}_fastqc_opd -f $trimfile"
    fastqc -o ${i}_fastqc_opd -f $trimfile 1>> process.log 2>&1
    #fastqc -o EO33_single_fastqc_opd -f fastq EO33_trim.fastq
    
    #bowtie
    echo "  Bowtie: Aligning $cleanfile to the genome using command:"
    echo "    bowtie -q --nomaqround --sam -m 4 -n 2 -e 70 -l 28 --best -p 8 --chunkmbs 1000 --seed=123 /proj/dllab/Erin/ce10/from_ucsc/seq/genome_bwa/ce10 $cleanfile $samfile"
    bowtie -q --nomaqround --sam -m 4 -n 2 -e 70 -l 28 --best -p 8 --chunkmbs 1000 --seed=123 /proj/dllab/Erin/ce10/from_ucsc/seq/genome_bwa/ce10 $cleanfile $samfile 1>> bt.log 2>&1
    #bsub -q week -n 8 -R span"[hosts=1]" -x -o bt.log bowtie -q --nomaqround --sam -m 4 -n 2 -e 70 -l 28 --best --chunkmbs 1000 -p 8 --seed=123 /proj/dllab/Erin/ce10/from_ucsc/seq/genome_bwa/ce10 EO23_clean.fastq EO23_output.sam
    #$cleanfile $samfile
    
    #samtools compress & sort
    echo "  Samtools:  Compressing $samfile into $bamfile using command:"
    echo "    samtools view -bS -o $bamfile $samfile"
    samtools view -bS -o $bamfile $samfile 1>> process.log 2>&1
    #bsub -q week -o log-samtools -n 8 -R span"[hosts=1]" -x samtools view -bS -o EO23.bam EO23_output.sam

    #samtools sort
    echo "  Samtools:  Sorting $bamfile into $bam_sort_file using command:"
    echo "    samtools sort $bamfile $bam_sorted"
    samtools sort $bamfile $bam_sorted 1>> process.log 2>&1
    #bsub -q week -o log-sort -n 8 samtools sort $bamfile $bam_sorted
    #bsub -q week -o log-sort -n 8 R span"[hosts=1]" -x samtools sort EO23.bam EO23_sorted
    
    #bedtools from .sam file
    echo "   Bedtools:  converting $bam_sort_file to $bed_file using command:"
    echo "     bedtools bamTobed -i $bam_sort_file > $bed_file"
    bedtools bamTobed -i $bam_sort_file > $bed_file
    #bedtools wig file
    echo "  R: writing a .R file to convert .bam --> 
    
    printf "\n"
    
#    bsub -q week -n 8 -o closest_EO33.LOG "bedtools closest -a EO37_overIgG_peaks.bed -b /proj/dllab/Erin/ce10/from_ucsc/annotation_files/WS220_refseq_ucsc_120323.gtf > EO37_closest.bed"
done


################
Notes:

#branching to output and log file
ls -l 2>&1 | tee file.txt

#zinba commands
library(zinba)
basealigncount(
  inputfile="EO37_overIgG.bed",
  outputfile="EO37.wig",
  extension=250,
  filetype="bed",
  twoBitFile="/proj/dllab/Erin/ce10/from_ucsc/seq/ce10.2bit"
)

library(zinba)
basealigncount(
  inputfile="EO27_sorted.bed",
  outputfile="EO27.wig",
  extension=250,
  filetype="bed",
  twoBitFile="/proj/dllab/Erin/ce10/from_ucsc/seq/ce10.2bit"
)

#Java outputs
/proj/dllab/Erin/executables/java-genomics-toolkit/toolRunner.sh wigmath.Subtract -f -m 
EO37.wig -o EO37_over_27 -s mEO27.wig>


#! /bin/bash/

logfilename="logfile1.log"

echo "$ wc EO38_mac3s.log..." > $logfilename
wc EO38_mac3s.log >> $logfilename 2>&1
echo "$ wc bt.log..." >> $logfilename
wc bt.log >> $logfilename 2>&1