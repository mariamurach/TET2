# TET2 manuscript code 

## QC and adapter/quality trimming

1. Fastqc - script like `qc.sh` but for untrimmed files
2. Use `fastp.sh` to remove adapters and filter based on quality
3. `qc.sh` - check quality, if adapters were removed

## Mapping

4. `star_index.sh` - create an index file for star
5. `map_array.sh` - map reads to reference
6. `samtools stats` qc of bam files to check the quality of mapping

## Retriving counts

7. `featureCounts_all.sh` - get counts of reads mapped to reference genes

## Differential analysis

8. `analysis/differential_expression.R` -  script for differential analysis using DEseq2 with firther pathways analysis. 

---
## Variables that need to be exported befor running the code: 
export curdir=/path/to/this/repo
export raw=$curdir/raw

## Notes

* Downstream bcr analysis is saved in https://github.com/mariamurach/bcr_R and uses output files from bcr_analysis.R

* Paths in the code files were anonymized and might require editing to fit users file stracture
* This pipeline uses code https://github.com/liulab-dfci/RIMA_pipeline to retrive CDR3 information and to concatenate bcr files from TRUST4 output. 
