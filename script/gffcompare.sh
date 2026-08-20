gffread Isoquant/OUT/OUT.transcript_models.gtf \
  -g Isoquant/GRCh38.primary_assembly.genome.fa \
  -w Isoquant/OUT.transcript_models.fa



# extracting the sequence for the novel transcripts

seqkit grep -r -f novel_transcripts_edger.txt Isoquant/OUT.transcript_models.fa -o novel_isoforms.fa


# running gffcompare

gffcompare -r Isoquant/GRCh38.gtf -o gffcmp_out Isoquant/OUT/OUT.transcript_models.gtf

# extract information related to novel transcripts only

grep -F -f novel9.txt Isoquant/OUT/gffcmp_out.OUT.transcript_models.gtf.tmap > novel9_classification.tsv


##explain what does the following commands do

gffread Isoquant/OUT/OUT.transcript_models.gtf \
  -g Isoquant/GRCh38.primary_assembly.genome.fa \
  -w Isoquant/OUT.transcript_models.fa


seqkit grep -f  novel_ids.txt OUT.transcript_models.fa -o novel_isoforms.fa



Q--> do the long noncoding RNA share sequence with the genes they regulate?
Q--> shall I use intersect or genomicRanges in R to compare the coordinates between the genes and my transcripts?
Q--> how can transdecoder give information about frameshits, premature stop codon, and novel start codon?
