#!/bin/bash

# Script to extract unique transcript_id from GTF files
# Usage: ./extract_transcript_ids.sh

# Set the directory containing GTF files
GTF_DIR="/home/vmaracaja/shared-projects/SUZ_sr_2025/exons"

# Create log file with timestamp
LOG_FILE="$GTF_DIR/transcript_extraction_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "$1" | tee -a "$LOG_FILE"
}

# Change to the GTF directory
cd "$GTF_DIR" || exit 1

# Start logging
log_message "=== Transcript ID Extraction Log ==="
log_message "Started at: $(date)"
log_message "Working directory: $GTF_DIR"
log_message ""

# Process each GTF file
for gtf_file in *.gtf; do
    if [ -f "$gtf_file" ]; then
        log_message "Processing: $gtf_file"
        
        # Extract basename without extension
        basename="${gtf_file%.gtf}"
        
        # Output file name
        output_file="${basename}_IDs.txt"
        
        # Extract unique transcript_id values
        # Look for transcript_id attribute in GTF files
        grep -o 'transcript_id "[^"]*"' "$gtf_file" | \
        sed 's/transcript_id "//g' | \
        sed 's/"//g' | \
        sort -u > "$output_file"
        
        # Count and display results
        count=$(wc -l < "$output_file")
        log_message "  -> Found $count unique transcript IDs in $output_file"
        log_message ""
    fi
done

log_message "All GTF files processed successfully!"
log_message "Completed at: $(date)"
log_message "Log saved to: $LOG_FILE"