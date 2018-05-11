## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/index.html)
- [bismark v0.19.0](https://github.com/FelixKrueger/Bismark)
- [R v3.3.2](https://www.r-project.org/)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)



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
