

for sample in FASTQ/*.fastq.gz
do
sample_name=$(basename "$sample" ".fastq.gz")
mkdir -p minimap_align/
minimap2 -t 12 -ax splice -uf -k14 reference/GRCh38.primary_assembly.genome.fa.gz "$sample" | samtools sort \
    | samtools view -@ 12 -bh > minimap_align/"$sample_name".bam
    samtools index minimap_align/"$sample_name".bam
done


> 

minimap_align/"$sample_name".sam
done





minimap2 -ax splice -uf -k14 reference/GRCh38.primary_assembly.genome.fa.gz FASTQ/PAU84062_pass_FAST_barcode01.fastq.gz > aln.sam  # noisy Nanopore Direct RNA-seq

