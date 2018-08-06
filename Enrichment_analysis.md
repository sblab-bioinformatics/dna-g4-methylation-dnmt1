##Enrichment analysis

ENCODE DHS and ChIP-seq data sets were normalised to sequencing depth of 1 (RPGC, Reads Per Genomic Content). Sequencing depth is defined as: (total number of mapped reads * fragment length) / effective genome size. The effective genome size was set to be 3,209,286,105 and enrichment values for DHSs and BG4 peaks over CGIs and their flanking sequences were visualised in R using ggplot2 library. Enrichment values for DNMT1 over CGIs and their flanks were visualised with DeepTools53.
