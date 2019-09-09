#!/bin/bash


if [ $# -ne 4 ]
then
    echo 'Usage: ./cellranger atac.sh [id] [fastqpath] [sampleID] \
    [t]'
    exit 65

fi     

id=$1
fastqpath=$2
sampleID=$3
t=$4


cellranger-atac count --id=$id --fastqs=$fastqpath \
                --reference=/share/data0/refdata-cellranger-atac-GRCh38-1.1.0 \
                --sample=$sampleID \
                --localcores=$t



              