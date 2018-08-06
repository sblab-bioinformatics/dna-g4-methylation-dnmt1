
## DNA G-quadruplex structures mould the DNA methylome

This repository contains data access and computational analysis for the methods developed in our manuscript *currently under revision*.

### Data

All the sequencing data have been deposited in the NCBI GEO database under accession number [GSE107690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107690).

### Code

- [**Whole genome bisulfite sequencing in HaCaT cells**](wgbs_hacat.md)
  - [Read processing, alignment and de-duplication](wgbs_hacat.md#read-processing-alignment-and-de-duplication)
  - [Methylation extraction, aggregation and filter by coverage](wgbs_hacat.md#methylation-extraction-aggregation-and-filter-by-coverage)
  - [Analysis of methylation changes in BG4 ChIP-Seq peaks, CGIs and promoters before and after Entinostat treatment](wgbs_hacat.md#analysis-of-methylation-changes-in-bg4-chip-seq-peaks-cgis-and-promoters-before-and-after-entinostat-treatment)
  
- **BG4 ChIP-seq in K562 cells**, following our [previous work](https://github.com/sblab-bioinformatics/dna-secondary-struct-chrom-lands/blob/master/Methods.md)

### Tools 

|Tool           | version  | URL (If available)                                            |
| ------------- |:--------:| --------------------------------------------------------------|
| cutadapt      | 1.15     |NA                                                             |
| BWA           | 0.7.15   |NA                                                             |
| Picard        | 2.8.3    |http://broadinstitute.github.io/picard                         |
| MACS          | 2.1.1    |NA                                                             |
| Bedtools      | 2.26.0   |http://bedtools.readthedocs.io/en/latest/content/overview.html |
| Deeptools     | 0.7.15   |NA                                                             |
| Bismark       | 0.19.0   |NA                                                             |

