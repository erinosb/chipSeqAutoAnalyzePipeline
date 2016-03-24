#! /bin/bash

for i in $@
do
    root=("${i%%.bed}")
    bedtools bedtobam -i ${i} -g /proj/dllab/Erin/ce10/from_ucsc/seq/chr_length_ce10.txt > ${root}.bam
    samtools sort ${root}.bam ${root}_sorted
    samtools index ${root}_sorted.bam
done