python sqanti3_qc.py \
  --isoforms ../novel_26_transcripts.gtf \
  --refGTF ../Isoquant/gencode.v47.annotation.gtf \
  --refFasta ../reference/GRCh38.primary_assembly.genome.fa \
  --SR_bam ../Star_aligned_sub/ \
  -c "../Star_aligned_sub/*/SJ.out.tab" \
  --fl_count ../fl_counts.tsv \
  --include_ORF \
  --aligner_choice minimap2 \
  -o sqanti3_novel26_subSR \
  -d sqanti3_out_subSR/ \
  -t 8 \
  --report  
    


python sqanti3_filter.py rules \
  --sqanti_class sqanti3_out_subSR/sqanti3_novel26_subSR_classification.txt \
  --filter_isoforms sqanti3_out_subSR/sqanti3_novel26_subSR_corrected.fasta \
  -o sqanti3_novel26_subSR_filtered \
  -d sqanti3_out_subSR/

