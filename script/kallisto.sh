#Kallisto


#index

 kallisto index -i index.idx reference/gencode.v47.transcripts.fa.gz


For sample in FASTQ/*.fastq-gz; do
    sample_name=$(basename "$sample" ".fastq.gz")
    echo "Working on $sample..."
	mkdir -p kallisto_out/"$sample_name"/
 kallisto bus -t 12 --long --threshold 0.8 -x bulk -i reference/ref_index \
  -o kallisto_out/"$sample_name"/ -G reference/gencode.v47.annotation.gtf "$sample"
done


for sample in kallisto_out/PAU*/;do
 bustools sort -t 12 "$sample"/output.bus \
 -o "$sample"/sorted.bus
done


for sample in kallisto_out/PAU*/;do
 bustools count "$sample"/sorted.bus \
 -t "$sample"/transcripts.txt \
 -e "$sample"/matrix.ec \
 -g reference/gencode.v47.annotation.gtf \
 -o "$sample"/count --cm -m
done
 

for sample in kallisto_out/PAU*/;do
kallisto quant-tcc -t 12 \
	--long -P ONT \
	 -f "$sample"/flens.txt \
	-i reference/ref_index \
	-e "$sample"/count.ec.txt \
	-G reference/gencode.v47.annotation.gtf \
	-o "$sample"/quant-tcc \
	--matrix-to-files \
	"$sample"/count.mtx
done

