#=======================================================================
# mapping rate results
#=======================================================================

import sys, os, glob

# ----------------------------------------------------------------------
# Set paths for log files and output.
# ----------------------------------------------------------------------
# 1. Directory path where the STAR log.final.out files are located.
LOG_DIRECTORY = "/path/to/STAR/log/files/"

# 2. Directory path to save the resulting CSV file (e.g., SW1783_mapping_rate.txt).
OUTPUT_DIRECTORY = "/path/to/output/results/"

# ----------------------------------------------------------------------
# 3. Process files and extract mapping statistics.
# ----------------------------------------------------------------------

# Define log file pattern and create file list
LOG_FILE_PATTERN = os.path.join(LOG_DIRECTORY, "*.final.out")
file_list = glob.glob(LOG_FILE_PATTERN)
file_list.sort() # Sort by file name

# Set output file path
OUTPUT_FILE_NAME = "sample_mapping_rate.txt"
output_file = os.path.join(OUTPUT_DIRECTORY, OUTPUT_FILE_NAME)

with open(output_file, 'w') as fw:
    # Write header
    header = "sample,Number of input reads,Uniquely mapped reads number,Uniquely mapped reads %\n"
    fw.write(header)
    
    for file_path in file_list:
        try:
            # Extract sample name from file path (e.g., sample.STAR.final.out -> sample)
            file_name = os.path.basename(file_path)
            samplename = file_name.split('.STAR')[0]
            
            with open(file_path, 'r') as f:
                # Read all lines from the file
                all_lines = list(f)
                
                # Extract values from fixed line indices in STAR log.final.out 
                # Values are split by '|' and tab (\t)
                
                # Line 5 (index 5): Number of input reads
                input_reads = all_lines[5].split('|\t')[1].strip()
                # Line 8 (index 8): Uniquely mapped reads number
                unique_reads = all_lines[8].split('|\t')[1].strip()
                # Line 9 (index 9): Uniquely mapped reads %
                mapping_rates = all_lines[9].split('|\t')[1].strip()
                
                # Use f-string to write the CSV record
                record = f"{samplename},{input_reads},{unique_reads},{mapping_rates}\n"
                fw.write(record)
        finally:
            pass
