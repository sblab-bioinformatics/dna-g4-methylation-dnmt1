## Software requirements

- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/index.html)
- [bismark v0.19.0](https://github.com/FelixKrueger/Bismark)


## Read processing, alignment and de-duplication

### Adaptor trimming and base quality filtering

The folder `fastq/` contains the raw paired-end fastq files (R1 and R2) from the whole genome bisulfite sequencing in the HaCaT cell line.

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
