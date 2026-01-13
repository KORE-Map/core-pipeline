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
