# core-pipeline
This repository contains the core analysis pipeline code consistently used across all yearly versions of the project.

## License

This project is licensed under the MIT License.

**Note:**  
This code is freely available for academic and research use.  
Commercial use may require a separate license agreement.  
Please contact the maintainers for inquiries regarding commercial usage.

## RNA_processing_code.sh
Main RNA-seq processing pipeline - trimming (Trim galore), alignment (STAR), quantification (RSEM)

## Sequencing Mapping Results.py
Parses STAR Log.final.out files to summarize mapping statistics (input reads, uniquely mapped reads, mapping rate)

## RNAseq_survey_example.R
Code for extracting and summarizing statistics from preprocessed RNA-seq data

### sample_info_example.xlsx
Sample organization example for running RNAseq_survey_example.R

## data_analysis_code.R
Performs DEG identification, Gene set enrichment analysis (GSEA), and CMap-based comparative transcriptomic analysis between RNA-seq data and public Connectivity Map (CMap) profiles.

### RNAseq_results_example.Rdata
### 20230719_Repurposing_Hub_export.txt
Sample organization example for running data_analysis_code.R
