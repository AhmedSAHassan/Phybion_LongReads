#!/bin/bash


# Long reads alignment and Quantification (working)

flair align -g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa -r reads/*.fastq.gz

flair correct -q flair.aligned.bed -f references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf -g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa

flair collapse -g references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa -q flair_all_corrected.bed -r reads/*.fastq.gz --gtf references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf

flair quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa







#creating minimap2 index

minimap2 -d ref.mmi reference/GRCh38.primary_assembly.genome.fa.gz  


# align with flair

for sample in FASTQ/*.fastq.gz; do
    sample_name=$(basename "$sample" ".fastq.gz")
    echo "Working on $sample..."
	mkdir -p flair_aligned/"$sample_name"/
    flair align \
        --mm_index ref.mmi \
        -r "$sample" \
        --threads 12 \
        --output flair_aligned/"$sample_name"/"$sample_name"
done


