## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/index.html)
- [bismark v0.19.0](https://github.com/FelixKrueger/Bismark)
- [R v3.3.2](https://www.r-project.org/)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [GenomicFeatures v1.26.4](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
  - [ggplot2 v2.2.1](https://ggplot2.tidyverse.org/)
- [python v2.7.12](https://www.python.org/)
- [fastaRegexFinder.py](https://github.com/dariober/bioinformatics-cafe/blob/master/fastaRegexFinder.py)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [pigz v2.3.4](https://zlib.net/pigz/)
- [bedGraphToBigWig v4](https://www.encodeproject.org/software/bedgraphtobigwig/)



## Sequencing files ids

| untreated     | 10 uM entinostat treated |
|:-------------:|:------------------------:|
| SLX-14400_S4  | SLX-14400_S3             |
| SLX-14400_S2  | SLX-14400_S1             |
| SLX-13601_S2  | SLX-13601_S1             |
| SLX-13610_S2  | SLX-13610_S1             |



## Read processing, alignment and de-duplication

### Adaptor trimming and base quality filtering

The folder `fastq/` contains the raw paired-end fastq files (R1 and R2) from the whole genome bisulfite sequencing in the HaCaT cell line (see sequencing files ids).

`cutdapt` was used with the following options:

```bash
cd fastq
mkdir ../fastq_trimmed

for fq1 in *_R1_001.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%_R1_001.fastq.gz}
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 10 -O 1 -u 6 -u -1 -U 6 -U -1 -o ../fastq_trimmed/$fq1 -p ../fastq_trimmed/$fq2 $fq1 $fq2 > ../fastq_trimmed/$bname.txt
done
```


### Alignment

The reference genome was prepared using `bismark_genome_preparation`.

The trimmed paired-end fastq files were first aligned to the reference genome using `bismark` with the following options:

```bash
cd fastq_trimmed
mkdir ../bismark

ref=../reference/

for fq1 in *_R1_001.fastq.gz
do
  bname=${fq1%_R1_001.fastq.gz}
  fq2=${fq1/_R1_/_R2_}
  bismark --non_directional --unmapped -o ../bismark -p 8 $ref -1 $fq1 -2 $fq2
done
```

The resulting unmapped paired-end fastq files were re-aligned in single-end mode to increase mapping efficiency:

```bash
cd bismark
mkdir ../bismark_se

ref=../reference

for fq in *.fq.gz
do
  bname=${fq%.fq.gz}
  bismark --non_directional --unmapped -p 8 --gzip -o ../bismark_se --temp_dir ~/tmp $ref $fq
done
```


### De-duplication

The paired-end alignments were deduplicated using `deduplicate_bismark` as follows:

```bash
cd bismark

for id in SLX-14400_S4 SLX-14400_S3 SLX-14400_S2 SLX-14400_S1 SLX-13601_S2 SLX-13601_S1 SLX-13610_S2 SLX-13610_S1 
do
  bams=`echo ${id}_L00[1-4]_R1_001_bismark_bt2_pe.bam`
  deduplicate_bismark -p --output_dir . --bam --multiple $bams
done
```

The single-end alignments were also deduplicated:

```bash
cd bismark_se

for id in SLX-14400_S4 SLX-14400_S3 SLX-14400_S2 SLX-14400_S1 SLX-13601_S2 SLX-13601_S1 SLX-13610_S2 SLX-13610_S1
do
  bams=`echo ${id}_L00[1-4]_R[1-2]_001.fastq.gz_unmapped_reads_[1-2]_bismark_bt2.bam`
  deduplicate_bismark -s --output_dir . --bam --multiple $bams
done
```



## Methylation extraction, aggregation and filter by coverage

### Methylation extraction

`bismark_methylation_extractor` was used to extract methylation from deduplicated paired-end alignments as follows:

```bash
cd bismark
mkdir ../methylation

ref=../reference/

for bam in *.deduplicated.bam
do
  bname=${bam%_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bam}
  bismark_methylation_extractor -p --comprehensive -o ../methylation --gzip --parallel 8 --bedGraph --cytosine_report --genome_folder $ref $bam
done
```

Methylation was also extracted from deduplicated single-end alignments like:

```bash
cd bismark_se
mkdir ../methylation_se

ref=../reference/

for bam in *.deduplicated.bam
do
  bname=${bam%_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bam}
  bismark_methylation_extractor -s --comprehensive -o ../methylation_se --gzip --parallel 8 --bedGraph --cytosine_report --genome_folder $ref $bam
done
```


### Aggregation

#### Combining paired-end and single-end methylation counts

This was performed in the R programming lagnguage as follows:

```r
library(data.table)

#############
# untreated #
#############

first_pe_data <- fread('zcat methylation/SLX-14400_S4_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
first_se_data <- fread('zcat methylation_se/SLX-14400_S4_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
second_pe_data <- fread('zcat methylation/SLX-14400_S2_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
second_se_data <- fread('zcat methylation_se/SLX-14400_S2_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
third_pe_data <- fread('zcat methylation/SLX-13601_S2_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
third_se_data <- fread('zcat methylation_se/SLX-13601_S2_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
fourth_pe_data <- fread('zcat methylation/SLX-13610_S2_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
fourth_se_data <- fread('zcat methylation_se/SLX-13610_S2_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')

untreated_data <- rbind(first_pe_data, first_se_data, second_pe_data, second_se_data, third_pe_data, third_se_data, fourth_pe_data, fourth_se_data)[order(V1, V2)]
setnames(untreated_data, names(untreated_data), c('chr', 'start', 'end', 'pct_met', 'cnt_met', 'cnt_unmet'))

untreated_data_aggregate <- untreated_data[, .(cnt_met = sum(cnt_met), cnt_unmet = sum(cnt_unmet)), by = .(chr, start, end)]
nrow(untreated_data_aggregate) # 51825458

untreated_data_aggregate[, start := start - 1]
untreated_data_aggregate[, end := end - 1]

gz <- gzfile("merge_all_2/untreated.bismark.cov.gz", "w")
write.table(untreated_data_aggregate, gz, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
close(gz)



############################
# 10 uM entinostat treated #
############################

first_pe_data <- fread('zcat methylation/SLX-14400_S3_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
first_se_data <- fread('zcat methylation_se/SLX-14400_S3_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
second_pe_data <- fread('zcat methylation/SLX-14400_S1_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
second_se_data <- fread('zcat methylation_se/SLX-14400_S1_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
third_pe_data <- fread('zcat methylation/SLX-13601_S1_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
third_se_data <- fread('zcat methylation_se/SLX-13601_S1_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')
fourth_pe_data <- fread('zcat methylation/SLX-13610_S1_L001_R1_001_bismark_bt2_pe.multiple.deduplicated.bismark.cov.gz')
fourth_se_data <- fread('zcat methylation_se/SLX-13610_S1_L001_R1_001.fastq.gz_unmapped_reads_1_bismark_bt2.multiple.deduplicated.bismark.cov.gz')

treated_data <- rbind(first_pe_data, first_se_data, second_pe_data, second_se_data, third_pe_data, third_se_data, fourth_pe_data, fourth_se_data)[order(V1, V2)]
setnames(treated_data, names(treated_data), c('chr', 'start', 'end', 'pct_met', 'cnt_met', 'cnt_unmet'))

treated_data_aggregate <- treated_data[, .(cnt_met = sum(cnt_met), cnt_unmet = sum(cnt_unmet)), by = .(chr, start, end)]
nrow(treated_data_aggregate) # 51720691

treated_data_aggregate[, start := start - 1]
treated_data_aggregate[, end := end - 1]

gz <- gzfile("merge_all_2/treated.bismark.cov.gz", "w")
write.table(treated_data_aggregate, gz, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
close(gz)
```

#### Collapse counts to CpG sites

Obtain all CpG sites in hg19:

```bash
mkdir cpg_sites
cd cpg_sites

ref=../reference/genome.fa
fastaRegexFinder.py -f $ref -r CG --noreverse | cut -f1-3 | grep "^chr" | bedtools sort -i | awk -v OFS="\t" '$3 = $3 - 1' | pigz > hg19.allCpG.bed.gz

zcat hg19.allCpG.bed.gz | wc -l # 28217448 CpGs, 56434896 Cs and Gs in CpGs
```

Intersect with `treated.bismark.cov.gz` and `untreated.bismark.cov.gz`:

```bash
cd ../merge_all_2

bedtools intersect -a <(zcat ../cpg_sites/hg19.allCpG.bed.gz) -b <(zcat treated.bismark.cov.gz | tail -n +2) -sorted -wa -wb | \
bedtools merge -i - -c 7,8 -o sum,sum | \
pigz > treated.collapse.bed.gz

bedtools intersect -a <(zcat ../cpg_sites/hg19.allCpG.bed.gz) -b <(zcat untreated.bismark.cov.gz | tail -n +2) -sorted -wa -wb | \
bedtools merge -i - -c 7,8 -o sum,sum | \
pigz > untreated.collapse.bed.gz

zcat treated.collapse.bed.gz | wc -l # 1x, 26686940 (95%)
zcat treated.collapse.bed.gz | awk '$4+$5 > 4' | wc -l # 5x, 22546472 (80%)
zcat treated.collapse.bed.gz | awk '$4+$5 > 9' | wc -l # 10x, 16534935 (59%)

zcat untreated.collapse.bed.gz | wc -l # 1x, 26692320 (95%)
zcat untreated.collapse.bed.gz | awk '$4+$5 > 4' | wc -l # 5x, 22670454 (80%)
zcat untreated.collapse.bed.gz | awk '$4+$5 > 9' | wc -l # 10x, 17044566 (60%)
```

Continue with 5x. Intersect `treated.collapse.bed.gz` and `untreated.collapse.bed.gz` at 5x coverage of individual CG sites to generate a consensus file containing CG sites in common to both conditions.

```bash
bedtools intersect -a <(zcat treated.collapse.bed.gz | awk -v OFS="\t" '$4+$5 > 4') -b <(zcat untreated.collapse.bed.gz | awk -v OFS="\t" '$4+$5 > 4') -sorted -wo | cut -f1-5,9-10 | pigz > treated.untreated.collapse.bed.gz
zcat treated.untreated.collapse.bed.gz | wc -l # 5x, 21106307 (75%) - after intersecting treated and untreated

# the format for each CG site is:
# <chr>
# <start>
# <end>
# <count_met_treated>
# <count_unmet_treated>
# <count_met_untreated>
# <count_unmet_untreated>
```

Convert `treated.collapse.bed.gz`, `untreated.collapse.bed.gz` and `treated.untreated.collapse.bed.gz` into BigWig format using `bedGraphToBigWig`:

```bash
#calculate chrom sizes as defined in bedGraphToBigWig help
cat ../reference/genome.fa.fai | cut -f1-2 | grep "chr" > ../reference/genome.sizes
ref=../genome.sizes

#treated.collapse.bed.gz
zcat treated.collapse.bed.gz | awk -v OFS="\t" '{print $1, $2, $3+1, 100*$4/($4+$5)}' > treated.collapse.bedgraph
bedGraphToBigWig treated.collapse.bedgraph $ref treated.collapse.bw
chmod 755 treated.collapse.bw
rm treated.collapse.bedgraph

#untreated.collapse.bed.gz
zcat untreated.collapse.bed.gz | awk -v OFS="\t" '{print $1, $2, $3+1, 100*$4/($4+$5)}' > untreated.collapse.bedgraph
bedGraphToBigWig untreated.collapse.bedgraph $ref untreated.collapse.bw
chmod 755 untreated.collapse.bw
rm untreated.collapse.bedgraph

#treated.untreated.collapse.bed.gz "treated"
zcat treated.untreated.collapse.bed.gz | awk -v OFS="\t" '{print $1, $2, $3+1, 100*$4/($4+$5)}' > treated.collapse.5x.intersection.bedgraph
bedGraphToBigWig treated.collapse.5x.intersection.bedgraph $ref treated.collapse.5x.intersection.bw
chmod 755 treated.collapse.5x.intersection.bw
rm treated.collapse.5x.intersection.bedgraph

#treated.untreated.collapse.bed.gz "untreated"
zcat treated.untreated.collapse.bed.gz | awk -v OFS="\t" '{print $1, $2, $3+1, 100*$6/($6+$7)}' > untreated.collapse.5x.intersection.bedgraph
bedGraphToBigWig untreated.collapse.5x.intersection.bedgraph $ref untreated.collapse.5x.intersection.bw
chmod 755 untreated.collapse.5x.intersection.bw
rm untreated.collapse.5x.intersection.bedgraph
```

Remember, the above BigWig files are % methylation considering CpG sites as individual units.



## Analysis of methylation changes in BG4 ChIP-Seq peaks, CGIs and promoters before and after Entinostat treatment

The following files were obtained from Figure 2d in [PMID: 27618450](https://www.nature.com/articles/ng.3662):

- 2d_ATAC+_OQS+_G4-.csv
- 2d_ATAC+_OQS+_G4ChIP+.csv
- 2d_ATAC+_OQS+_G4ChIP++.csv

Click [here](https://media.nature.com/original/nature-assets/ng/journal/v48/n10/source_data/ng.3662-f2.xls) to download Figure 2d data.

```r
library(data.table)
library(GenomicFeatures)

# Load data from the paper
ATACp_OQSp_G4m <- data.table(read.csv("genes/2d_ATAC+_OQS+_G4-.csv"))
nrow(ATACp_OQSp_G4m) # 1734
ATACp_OQSp_G4m[, id := "ATACp_OQSp_G4m"]

ATACp_OQSp_G4ChIPp <- data.table(read.csv("genes/2d_ATAC+_OQS+_G4ChIP+.csv"))
nrow(ATACp_OQSp_G4ChIPp) # 3627
ATACp_OQSp_G4ChIPp[, id := "ATACp_OQSp_G4p"]

ATACp_OQSp_G4ChIPpp <- data.table(read.csv("genes/2d_ATAC+_OQS+_G4ChIP++.csv"))
nrow(ATACp_OQSp_G4ChIPpp) # 373
ATACp_OQSp_G4ChIPpp[, id := "ATACp_OQSp_G4pp"]

paper_tables <- rbindlist(list(ATACp_OQSp_G4m, ATACp_OQSp_G4ChIPp, ATACp_OQSp_G4ChIPpp))

# Prepare tss coordinates
txdb <- makeTxDbFromGFF("reference/genes.gtf", format="gtf")
genes_coords <- data.table(as.data.frame(genes(txdb)))
tss_1000_coords <- genes_coords[, .(chr = seqnames, start = ifelse(strand == "+", start - 1000, end - 1000), end = ifelse(strand == "+", start + 1000, end + 1000), gene_name = gene_id)]

# Merge
paper_tables_tss_1000_coords <- merge(paper_tables, tss_1000_coords, all.x = TRUE, by.x = "gene_id", by.y = "gene_name")
write.table(paper_tables_tss_1000_coords[!is.na(chr)][,.(chr, start, end, gene_id, logFC, FDR, id)], file = "genes/tss.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```

Obtain TSSs overlapping with CGIs (obtained from UCSC's cpgIslandExt.txt.gz table):

```bash
cd merge_all_2

zcat ../cpgi/cpgIslandExt.txt.gz | cut -f2-4 | tail -n +2 | bedtools sort -i > ../cpgi/cpgIslandExt.bed
wc -l ../cpgi/cpgIslandExt.bed # 28691

grep -v "_" ../cpgi/cpgIslandExt.bed > ../cpgi/cpgIslandExt_norandom.bed

cd ../genes

bedtools intersect \
-a tss.bed \
-b ../cpgi/cpgIslandExt_norandom.bed \
-wa -u > tss.cgi.bed

wc -l tss.cgi.bed # 5072
```

Obtain CpG sites located within the TSS regions obtained above:

```bash
cd merge_all_2

bedtools intersect \
-a <(zcat treated.untreated.collapse.bed.gz) \
-b <(bedtools sort -i ../genes/tss.cgi.bed) \
-sorted -wa -u | wc -l # 444526 CpG sites

bedtools intersect \
-a <(zcat treated.untreated.collapse.bed.gz) \
-b <(bedtools sort -i ../genes/tss.cgi.bed) \
-sorted -wb | cut -f8-10 | sort | uniq -c | wc -l # 5070 gene TSSs

# 444526 CpG sites overlap with 5070 gene TSSs.

bedtools intersect \
-a <(zcat treated.untreated.collapse.bed.gz) \
-b <(bedtools sort -i ../genes/tss.cgi.bed) \
-sorted -wao | cut -f1-14 | pigz > treated.untreated.collapse.tss.cgi.bed.gz

zcat treated.untreated.collapse.tss.cgi.bed.gz | wc -l # 21148404

# the format for each CpG site in treated.untreated.collapse.tss.bed.gz is:
# <chr>
# <start>
# <end>
# <count_met_treated>
# <count_unmet_treated>
# <count_met_untreated>
# <count_unmet_untreated>
# <chr overlapping tss> "." if there is no overlap
# <start overlapping tss> "-1" if there is no overlap
# <end overlapping tss> "-1" if there is no overlap
# <gene_id overlapping tss> "." if there is no overlap
# <logFC overlapping tss> "." if there is no overlap
# <FDR overlapping tss> "." if there is no overlap
# <group_id overlapping tss> "." if there is no overlap
```

Analysing `treated.untreated.collapse.tss.cgi.bed.gz` (as generated above) in R,

```r
library(data.table)
library(ggplot2)

# Load data
data <- fread("zcat merge_all_2/treated.untreated.collapse.tss.cgi.bed.gz")
colnames(data) <- c("chr_cpg", "start_cpg", "end_cpg", "cnt_met_treated", "cnt_unmet_treated", "cnt_met_untreated", "cnt_unmet_untreated", "chr_tss", "start_tss", "end_tss", "gene_tss", "logFC_tss", "FDR_tss", "id_tss")

# Calculate pct_met for each CpG and diff
data[, pct_met_treated := 100*cnt_met_treated/(cnt_met_treated+cnt_unmet_treated)]
data[, pct_met_untreated := 100*cnt_met_untreated/(cnt_met_untreated+cnt_unmet_untreated)]
data[, diff := 100*cnt_met_treated/(cnt_met_treated+cnt_unmet_treated) - 100*cnt_met_untreated/(cnt_met_untreated+cnt_unmet_untreated)]
data[, logFC_tss := as.numeric(data$logFC_tss)]

# Collapse data to TSS taking the sum of methylation counts, average of pct_met_treated, pct_met_untreated and diff across CpG sites, counting how many CpG sites we have for each intersection, then calculated significantly methylated TSSs with a Fisher's test on the counts
data_tss <- data[chr_tss != "."][, .(cnt_met_treated = sum(cnt_met_treated), cnt_unmet_treated = sum(cnt_unmet_treated), cnt_met_untreated = sum(cnt_met_untreated), cnt_unmet_untreated = sum(cnt_unmet_untreated), pct_met_treated = mean(pct_met_treated), pct_met_untreated = mean(pct_met_untreated), diff = mean(diff), .N), by = .(chr_tss, start_tss, end_tss, gene_tss, logFC_tss, FDR_tss, id_tss)]
prior <- 0.001
data_tss[, fc := ((cnt_met_treated/(cnt_met_treated+cnt_unmet_treated))+prior)/((cnt_met_untreated/(cnt_met_untreated+cnt_unmet_untreated))+prior)]
data_tss[, pval_f := fisher.test(matrix(c(cnt_met_treated, cnt_unmet_treated, cnt_met_untreated, cnt_unmet_untreated), nrow = 2))$p.value, by = 1:nrow(data_tss)]
data_tss[, pval_f_adj := p.adjust(data_tss$pval_f, method = "BH")]

table(data_tss$id_tss)
# ATACp_OQSp_G4m  ATACp_OQSp_G4p ATACp_OQSp_G4pp
#           1504            3261             307

summary(data_tss[, c("logFC_tss", "pct_met_treated", "pct_met_untreated", "diff", "fc", "N")])
#   logFC_tss        pct_met_treated  pct_met_untreated      diff               fc                N         
# Min.   :-5.17330   Min.   : 0.000   Min.   : 0.000    Min.   :-7.9049   Min.   :0.07653   Min.   : 11.00  
# 1st Qu.:-0.44920   1st Qu.: 2.252   1st Qu.: 2.304    1st Qu.:-0.6700   1st Qu.:0.95113   1st Qu.: 76.00  
# Median :-0.07632   Median : 6.696   Median : 6.834    Median :-0.1450   Median :1.10881   Median : 93.00  
# Mean   : 0.03388   Mean   : 8.986   Mean   : 9.158    Mean   :-0.1721   Mean   :1.15289   Mean   : 95.94  
# 3rd Qu.: 0.36178   3rd Qu.:13.046   3rd Qu.:13.348    3rd Qu.: 0.3199   3rd Qu.:1.30202   3rd Qu.:113.00  
# Max.   : 7.47543   Max.   :76.948   Max.   :76.619    Max.   :15.3822   Max.   :6.27281   Max.   :226.00  

data_tss[id_tss == "ATACp_OQSp_G4m", code := "-\n+\n+"]
data_tss[id_tss == "ATACp_OQSp_G4p", code := "+\n+\n+"]
data_tss[id_tss == "ATACp_OQSp_G4pp", code := "++\n+\n+"]

# diff
gg <- ggplot(data_tss, aes(x = code, y = diff)) +
geom_hline(yintercept = 0, linetype = "dotted") +
geom_jitter(colour = 'gray', alpha = 0.35, width = 0.3, size = 0.2) +
geom_boxplot(outlier.shape=NA, alpha = 0) +
coord_cartesian(ylim = c(-3, 3)) +
ylab("% methylation difference (treated - untreated)") +
xlab("") +
theme_classic() +
theme(axis.title=element_text(size=16), axis.text.y=element_text(size=16), axis.text.x=element_text(size=11))

ggsave('figures/20180423_treated.untreated.collapse.tss.cgi.boxplot.diff.gray.png', width = 12, units= 'cm')

mean(data_tss[id_tss == "ATACp_OQSp_G4m"]$diff) # -0.1426609
mean(data_tss[id_tss == "ATACp_OQSp_G4p"]$diff) # -0.1665028
mean(data_tss[id_tss == "ATACp_OQSp_G4pp"]$diff) # -0.3763876

wilcox.test(data_tss[id_tss == "ATACp_OQSp_G4pp"]$diff, data_tss[id_tss == "ATACp_OQSp_G4m"]$diff, alternative = "less")$p.value # 0.0003326281
wilcox.test(data_tss[id_tss == "ATACp_OQSp_G4pp"]$diff, data_tss[id_tss == "ATACp_OQSp_G4p"]$diff, alternative = "less")$p.value # 1.301178e-05
wilcox.test(data_tss[id_tss == "ATACp_OQSp_G4p"]$diff, data_tss[id_tss == "ATACp_OQSp_G4m"]$diff, alternative = "less")$p.value # 0.5075517

kruskal.test(diff ~ as.factor(id_tss), data = data_tss)$p.value # 0.0002611352

pairwise.wilcox.test(data_tss$diff, data_tss$id_tss, p.adjust.method = "BH")
#                ATACp_OQSp_G4m ATACp_OQSp_G4p
#ATACp_OQSp_G4p  0.985          -             
#ATACp_OQSp_G4pp 0.001          7.8e-05       
```
