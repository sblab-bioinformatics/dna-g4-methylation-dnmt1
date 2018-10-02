
## DNA G-quadruplex structures mould the DNA methylome

This repository contains data access and computational analysis for the methods developed in our [published manuscript](https://www.nature.com/articles/s41594-018-0131-8).

### Data

- All in-house the sequencing data have been deposited in the NCBI GEO database under accession number [GSE107690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107690). 

- [**ENCODE**](https://www.encodeproject.org/)'s K562 cell line datasets used are as below: 

| Dataset         | Data Type                          | URL                                                    |
| --------------- |:----------------------------------:| -------------------------------------------------------|
| K562 DHSs       | DNase-seq                          | https://www.encodeproject.org/experiments/ENCSR000EPC/ |
| K562 DNMT1      | ChIP-seq                           | https://www.encodeproject.org/experiments/ENCSR987PBI/ |
| K562 WGBS       | Whole Genome Bisulfite Sequencing  | https://www.encodeproject.org/experiments/ENCSR765JPC/ |

- CGI were downloaded using the UCSC’s table browser. Also, as an alternative CpG Island definition, we used [CpGCluster](http://bioinfo2.ugr.es:8080/CpGislandDB/)

### Code

- [**BG4 ChIP-seq analysis**](ChIP-seq_Analysis.md)
- [**Analyses of methylation at promoters and CpG Islands**](Analyses_of_methylation_at_promoters_and_CpG_Islands.md)
- [**Enrichment analysis**](Enrichment_analysis.md)
- [**Monte Carlo Simulation**](Monte_Carlo_Simulation.md)
  - [Monte Carlo P-value calculation](Monte_Carlo_Simulation.md#Monte-Carlo-P-value-Calculation)
- [**Whole genome bisulfite sequencing in HaCaT cells**](wgbs_hacat.md)
  - [Read processing, alignment and de-duplication](wgbs_hacat.md#read-processing-alignment-and-de-duplication)
  - [Methylation extraction, aggregation and filter by coverage](wgbs_hacat.md#methylation-extraction-aggregation-and-filter-by-coverage)
  - [Analysis of methylation changes in BG4 ChIP-Seq peaks, CGIs and promoters before and after Entinostat treatment](wgbs_hacat.md#analysis-of-methylation-changes-in-bg4-chip-seq-peaks-cgis-and-promoters-before-and-after-entinostat-treatment)
  


### Tools 

|Tool           | version                         | URL                                                           |
| ------------- |:--------------------------------| --------------------------------------------------------------|
| cutadapt      | 1.15                            |http://cutadapt.readthedocs.io/en/stable/installation.html     |
| BWA           | 0.7.15                          |http://bio-bwa.sourceforge.net/                                |
| Picard        | 2.8.3                           |http://broadinstitute.github.io/picard                         |
| MACS          | 2.1.1                           |http://liulab.dfci.harvard.edu/MACS/                           |
| Bedtools      | 2.26.0                          |http://bedtools.readthedocs.io/en/latest/content/overview.html |
| Deeptools     | 0.7.15                          |https://deeptools.readthedocs.io/                              |
| Bismark       | 0.19.0                          |https://www.bioinformatics.babraham.ac.uk/projects/bismark/    |
| Liftover      | Accessed Aug 2017 to April 2018 |http://genome.ucsc.edu/cgi-bin/hgLiftOver                      |

