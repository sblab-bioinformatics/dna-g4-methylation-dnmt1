## Analyses of methylation at promoters and CpG Islands using ENCODE WGBS dataset

- Promoters were defined as 1 kb (+/−) from the transcription start sites of 31,239 hg19 transcripts. 
- Methylation levels at CpG sites with less than 5x coverage were discarded. 
- CGI were downloaded using the UCSC’s table browser and then ported to human genome release hg38 using the batch coordinate conversion (liftover) tool of the UCSC. 
- The alternative CGI sets were generated using [CpGCluster](http://bioinfo2.ugr.es:8080/CpGislandDB/). 
