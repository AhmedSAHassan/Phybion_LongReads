#!/bin/bash

#Long read quality control (working)

for sample in *.fastq.gz
do
sample_name=$(basename "$sample" ".fastq.gz")
mkdir -p nanoplot/"$sample_name"
NanoPlot \
--fastq_rich "$sample" \
--outdir nanoplot/"$sample_name"
done
