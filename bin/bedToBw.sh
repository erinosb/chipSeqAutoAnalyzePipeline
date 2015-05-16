#!/bin/bash

# Copyright (c) 2015, Erin Nishimura
#
# All rights reserved.
# Adapted from Kohta's 2012 script: bedToBw.sh 
#
#
#EXAMPLE:
#bash bedToBw.sh 7_EO40.bed 250 /proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt -n -bw
#
# April 20, 2015 -- Updated to accommodate a .bed file.

    printf "bedToBw:   START\n"
    
if [ $# -lt 5 ] # '-lt' stands for less than 
then
	echo -e "\t`basename $0` version 1"
	echo -e "\tby Erin Nishimura; all rights reserved"
  	echo -e "\n\tUsage: `basename $0` <file.bed> <extension> <Genome> -c/n -r/bg"
	echo -e "\n\t<hits.bed>\tBowtie output. Quality is not considered"
	echo -e "\t<extension>\tDesired final fragment size in bp."
	echo -e "\t<genome>\tSorted chromosome length file."
	echo -e "\t-c/n\t\tBase count (-c) or coverage-normalized base count (-n)."
	echo -e "\t-r/bg/bw\t\tOutput is read length only or bedgraph."
	echo -e "\n\tstdout not available\n"
	echo -e "\n\tRequires bedtools preinstalled\n"
	exit 
fi

# USAGE:
# $0 function name
# $1 hits.bed file in which each entry is a mapped read
# $2 Extension (final length)
# $3 Genome
# $4 Normalize -c/n
# $5 Output format: -r Print only read length; -bg bedGraph; -bw bigwig



# 1) Determine output file names
	F=$1
	seed="${F%%.bed}"	
	if [ $4 = "-c" ]; then	
		bgout="$seed"x"$2"c.bg
		bwout="$seed"x"$2"c.bw
	elif [ $4 = "-n" ]; then	
		bgout="$seed"x"$2"n.bg
		bwout="$seed"x"$2"n.bw
	else 
		echo -e "\tNormalization option (-n/c) required. Exit."
		exit
	fi
        
        #echo $seed
        #printf "\n"
	
        
# 2) Get read length. Get this from the wc of the input bed file ($1)
	rlength=$(head $1 -n 1 | awk '{print($3-$2)}') 
	
	if [ $5 = "-r" ]; then
		echo -e "bedToBw:   Here is the read length:" $rlength
		exit
	fi
        
        printf "bedToBw:   Read length is $rlength \n"
        
        
# 3) Get the genome length. Get this from the sum of the chromosome length file ($3)
        
        printf "bedToBw:   Get genome size\n"
	
        genome=$(cat $3 | awk '{total+=$2}END{print total}')

        echo -e "bedToBw:   Genome size is: $genome"

# 4) Calculate scale factor. If $4 is "-c", scale=1. If $4 is "-n", scale factor is GenomeSize / (final#BpExtendedReads * #ReadsInInputBedFile)
	if [ $4 = "-c" ]; then
		scale=1
	elif [ $4 = "-n" ]; then
                mappedreads=$(wc -l $1 | awk '{ print $1}')
                scale=$(wc -l $1 | awk '{printf "%.2f", '$genome'/('$2'*$1)}')
                echo -e "bedToBw:   The number of mapped reads is: ${mappedreads}."
	else
		echo -e "\tNormalization option (-n/c) required. Exit."
		exit
	fi
        
        echo -e "bedToBw:   Normalization factor is: $scale"

# 5) Sort using "bedtools sort"
# 6) Extend reads by final length minus the read length ($2 - $rlength)
# 7) Convert to a bedgraph using "genomeCoverageBed" (bedtools)
        final=$2
        n="$((final - rlength))"
        printf "bedToBw:   Bed file extend by $n bases.\n"
        
        printf "bedToBw:   Sorting\n"
	bedtools sort -i ${1} | bedtools slop -s -i stdin -g ${3} -l 0 -r ${n} | bedtools genomecov -bg -i stdin -g ${3} -scale ${scale} > ${bgout}


# 8) Convert bedgraph to a bigwig with bedGraphToBigWig:
        printf "bedToBw:   Converting to bigwig\n"
        
        /proj/.test/roach/FAIRE/bin/bedGraphToBigWig ${bgout} ${3} ${bwout}

# 9) cleanup ${bgout}
        
        rm ${bgout}
        
        printf "bedToBw:   END\n"
