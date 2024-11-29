
isoquant.py --reference reference/GRCh38.primary_assembly.genome.fa.gz \
  --genedb reference/gencode.v47.annotation.gtf \
  --bam minimap_align/PAU84062_pass_FAST_barcode01/PAU84062_pass_FAST_barcode01.bam minimap_align/PAU84062_pass_FAST_barcode02/PAU84062_pass_FAST_barcode02.bam minimap_align/PAU84062_pass_FAST_barcode03/PAU84062_pass_FAST_barcode03.bam minimap_align/PAU84062_pass_FAST_barcode04/PAU84062_pass_FAST_barcode04.bam minimap_align/PAU84062_pass_FAST_barcode05/PAU84062_pass_FAST_barcode05.bam minimap_align/PAU84062_pass_FAST_barcode06/PAU84062_pass_FAST_barcode06.bam\
  --data_type nanopore -o Isoquant/ --complete_genedb\
  --model_construction_strategy default_ont --report_novel_unspliced true


isoquant.py -d nanopore --bam ONT.cDNA_1.bam ONT.cDNA_2.bam ONT.cDNA_3.bam \
 --reference reference.fasta --genedb annotation.gtf --complete_genedb --output output_dir
 --predix ONT_3samples --labels A1 A2 A3

isoquant.py -d nanopore --yaml dataset.yaml  \
 --complete_genedb --genedb genes.gtf \
 --reference reference.fasta --output output_dir