#!/bin/bash

#R1----cell barcode and UMI 16+10bp --r1-length must >= 26
#R2----cDNA length
if [ $# -ne 6 ]
then
    echo "Usage: ./cellranger.sh [id] [fastqpath] [sampleID] \
    [t] [r1_length] [r2_length]"
    exit 65
fi

id=$1
fastqpath=$2
sampleID=$3
t=$4
r1_length=$5
r2_length=$6


cellranger count --id=$id --fastqs=$fastqpath \
	--transcriptome=/share/swap/refdata-cellranger-hg19-1.2.0/ \
	--sample=$sampleID --localcores=$t \
	--r1-length=$r1_length --r2-length=$r2_length

