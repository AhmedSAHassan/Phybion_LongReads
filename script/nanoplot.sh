#!/bin/bash

#Long read quality control (working)

for sample in *.fastq
do
sample_name=$(basename "$sample" ".fastq")
mkdir -p nanoplot/"$sample_name"
NanoPlot \
--fastq_rich "$sample" \
--outdir nanoplot/"$sample_name"
done
