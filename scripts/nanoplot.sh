

#Long read quality control

for sample in FASTQ/*.fastq.gz
do
sample_name=$(basename "$sample" ".fastq.gz")
mkdir -p nanoplot_SUB/"$sample_name"
NanoPlot \
--fastq "$sample" \
--outdir nanoplot_SUB/"$sample_name"
done
