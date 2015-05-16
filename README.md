# chipSeqAutoAnalyzePipeline

This is the basic initial pipeline we use for analyzing multiplexed ChIP-seq datasets

There are three steps to executign the analysis:
    chipSeqAnalyzeStep1.sh
    chipSeqAnalyzeStep2.sh


# chipSeqAnalyzeStep1.sh
  Copyright (c) 2015, Erin Osborne Nishimura

## PROGRAM
   autoAnalyzeChipseq_v9.sh - To automate the quality control & alignment of multiplexed ChIP-seq data


## USAGE
   Split and Align Mode:
       bash autoAnalyzeChipseq_v9.sh [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
   OR
   Split Only Mode:
       bash autoAnalyzeChipseq_v9.sh --splitOnly [options] --multi <inputFile.txt> --bar <barcodeIndexFile.txt>
   OR
   Align Only Mode:
       bash autoAnalyzeChipseq_v9.sh --alignOnly [options] <input1.fastq> input2.fastq


## MODES
   There are three modes for use... --splitNAlign mode (default), --splitOnly and --alignOnly mode
   --splitNAlign (default)      This mode uses as input a single homebrew multiplexed .txt file and
                               a single barcode index file. The pipeline then...
                                   0) splits the multiplexed file into single sample files by index
                                   1) Use fastx_trimmer to remove the barcode index from each sequencing read and remove low quality reads
                                   2) perform Tagdust using a solexa library of adapter and primer sequences.
                                   3) perform Fastqc to make a report of quality
                                   4) bowtie alignment
                                   5) compress .sam --> .bam using samtools view
                                   6) Sort .bam into \_sorted.bam using samtools sort
                                   7) Convert \_sorted.bam into .bed using bedtools
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
                                   6) Sort .bam into \_sorted.bam using samtools sort
                                   7) Convert \_sorted.bam into .bed using bedtools
                                   8) Convert .bed into .bw file using bedToBw.sh dependency script


## ARGUMENTS
   <inputFile.txt>           This is the file from the sequencing facility. It is a fastq file containing Illumina sequencing reads generated from multiplexed samples
   <barcodeIndexFile.txt>    This is a file listing the barcodes that are at the 5' end of each sequence string.  It should be in the format:

                   #barcode file for the multiplexed sequences generated 11/11/12
                   EO33	AGATGGT
                   EO34	CACGTCG
                   EO35	GATCTTG
                   EO36	TCAGGAC
                   EO37	ACAGTTG
   <input1.fastq>          This is an input .fastq file 


## OPTIONS

   --splitOnly                 Runs in splitOnly mode. Suppresses analysis.
   --alignOnly                 Runs in alignOnly mode. Suppresses splitting sequencefiles
   --qualityOff                Suppresses quality score reports
   --cleanOff                  Runs without the cleanup mode loop at the very end. The cleanup mode loop removes the \_clean.fastq, .sam, .bam, and .bed files.
                                 It retains the split.fastq.gz, \_trim.fastq.gz, \_sorted.bam.gz, and .bw files
   --trimOff                   Suppresses trimming six basepairs of multiplexing indices

   --extension <n>             Specify the length bp to extend reads in the .wig file. Default = 100
   -p <n>                      Runs bowtie in parallel mode. Values accepted are 1 - 8. Suggest 2 - 4. Make sure to match this number with bsub -n <n>. Default = 1

   --chrlength                 The bedToBw.sh program needs to know how long each chromsome is. Specify the location of the chromosome length file. Default =
                                       "/proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt"
   --bowtiepath                Bowtie1 requires the path location of where the .bwq files are contained. Default =
                                       "/proj/dllab/Erin/ce10/from_ucsc/seq/prev_versions_bowtie/genome_bwa/ce10"
   --primeradapter             Tagdust needs a .fastq file that contains a list of all the primer and adpater sequences used in library prep and sequencing. Tagdust
                                 remove these sequences. Default =
                                       "/proj/dllab/Erin/sequences/solexa-library-seqs.fasta"
   

## DEPENDENCIES
   Requires tagdust, fastqc, fastx-toolkit, bowtie1, bedtools, samtools, bedGraphToBigWig

   Developed with versions: TagDust 1.12; FastQC v0.11.3; bowtie/1.1.0; bedtools/2.22.1; samtools/0.1.19 


## AUTHOR
   Erin Osborne Nishimura


## DATE
   May 16, 2015

## BUGS/FUTURE EXPANSION
   -- check whether certain modules have been loaded;
   -- Auto load all required modules
   -- toggle between bowtie or bowtie2