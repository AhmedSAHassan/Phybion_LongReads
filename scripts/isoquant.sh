
isoquant.py --reference reference/GRCh38.primary_assembly.genome.fa.gz \
  --genedb reference/gencode.v47.annotation.gtf \
  --bam minimap_align/PAU84062_pass_FAST_barcode01/PAU84062_pass_FAST_barcode01.bam minimap_align/PAU84062_pass_FAST_barcode02/PAU84062_pass_FAST_barcode02.bam minimap_align/PAU84062_pass_FAST_barcode03/PAU84062_pass_FAST_barcode03.bam minimap_align/PAU84062_pass_FAST_barcode04/PAU84062_pass_FAST_barcode04.bam minimap_align/PAU84062_pass_FAST_barcode05/PAU84062_pass_FAST_barcode05.bam minimap_align/PAU84062_pass_FAST_barcode06/PAU84062_pass_FAST_barcode06.bam\
  --data_type nanopore -o Isoquant/ --complete_genedb\
  --model_construction_strategy default_ont --report_novel_unspliced true







