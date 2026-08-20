gffread Isoquant/OUT/OUT.transcript_models.gtf \
  -g Isoquant/GRCh38.primary_assembly.genome.fa \
  -w Isoquant/OUT.transcript_models.fa




# running gffcompare

gffcompare -r Isoquant/GRCh38.gtf -o gffcmp_out Isoquant/OUT/OUT.transcript_models.gtf




