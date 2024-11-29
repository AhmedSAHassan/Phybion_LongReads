for sample in minimap_align/PAU*; do
    echo "Working on $sample..."
sample_name=$(basename "$sample")
oarfish -j 16 -a "$sample"/"$sample_name".bam -o quants/"$sample_name" --filter-group no-filters --model-coverage
done
