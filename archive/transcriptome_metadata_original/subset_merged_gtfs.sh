#!/bin/bash

# Script to subset merged GTF files using transcript IDs
# Usage: ./subset_merged_gtfs.sh

# Set directories
MAIN_DIR="/home/vmaracaja/shared-projects/SUZ_sr_2025"
EXONS_DIR="$MAIN_DIR/exons"
TRANSCRIPTS_DIR="$MAIN_DIR/transcripts"

# Create transcripts directory if it doesn't exist
mkdir -p "$TRANSCRIPTS_DIR"

# Create log file with timestamp
LOG_FILE="$MAIN_DIR/gtf_subset_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log_message() {
    echo "$1" | tee -a "$LOG_FILE"
}

# Change to the main directory
cd "$MAIN_DIR" || exit 1

# Start logging
log_message "=== GTF Subsetting Log ==="
log_message "Started at: $(date)"
log_message "Working directory: $MAIN_DIR"
log_message "Output directory: $TRANSCRIPTS_DIR"
log_message ""

# Process each ID file in the exons directory
for id_file in "$EXONS_DIR"/*_IDs.txt; do
    if [ -f "$id_file" ]; then
        # Extract basename from ID file (remove path and _IDs.txt)
        id_basename=$(basename "$id_file" _IDs.txt)
        
        # Determine the corresponding merged GTF file
        merged_gtf=""
        output_gtf=""
        
        # Map ID files to merged GTF files
        case "$id_basename" in
            "CNE_24_HISAT2.annotated")
                merged_gtf="merged_CNE_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_HISAT2.annotated_tx.gtf"
                ;;
            "CNE_24_HISAT2.annotated.lncRNA")
                merged_gtf="merged_CNE_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_HISAT2.annotated.lncRNA_tx.gtf"
                ;;
            "CNE_24_HISAT2.annotated.otherRNAs")
                merged_gtf="merged_CNE_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_HISAT2.annotated.otherRNAs_tx.gtf"
                ;;
            "CNE_24_STAR.annotated.lncRNA")
                merged_gtf="merged_CNE_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_STAR.annotated.lncRNA_tx.gtf"
                ;;
            "CNE_24_STAR.annotated.otherRNAs")
                merged_gtf="merged_CNE_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_STAR.annotated.otherRNAs_tx.gtf"
                ;;
            "CNE_6_24_HISAT2.annotated")
                merged_gtf="merged_CNE_6_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_HISAT2.annotated_tx.gtf"
                ;;
            "CNE_6_24_HISAT2.annotated.lncRNA")
                merged_gtf="merged_CNE_6_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_HISAT2.annotated.lncRNA_tx.gtf"
                ;;
            "CNE_6_24_HISAT2.annotated.otherRNAs")
                merged_gtf="merged_CNE_6_24_HISAT2.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_HISAT2.annotated.otherRNAs_tx.gtf"
                ;;
            "CNE_6_24_STAR.annotated")
                merged_gtf="merged_CNE_6_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_STAR.annotated_tx.gtf"
                ;;
            "CNE_6_24_STAR.annotated.lncRNA")
                merged_gtf="merged_CNE_6_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_STAR.annotated.lncRNA_tx.gtf"
                ;;
            "CNE_6_24_STAR.annotated.otherRNAs")
                merged_gtf="merged_CNE_6_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_6_24_STAR.annotated.otherRNAs_tx.gtf"
                ;;
            "CNE24_STAR.annotated")
                merged_gtf="merged_CNE_24_STAR.annotated.gtf"
                output_gtf="$TRANSCRIPTS_DIR/CNE_24_STAR.annotated_tx.gtf"
                ;;
            *)
                log_message "Warning: No mapping found for $id_basename"
                continue
                ;;
        esac
        
        # Check if merged GTF file exists
        if [ ! -f "$merged_gtf" ]; then
            log_message "Warning: Merged GTF file $merged_gtf not found for $id_basename"
            continue
        fi
        
        log_message "Processing: $id_basename"
        log_message "  ID file: $id_file"
        log_message "  Merged GTF: $merged_gtf"
        log_message "  Output GTF: $output_gtf"
        
        # Count IDs to subset
        id_count=$(wc -l < "$id_file")
        log_message "  Subsetting with $id_count transcript IDs"
        
        # Create a temporary file with the pattern for grep
        temp_pattern_file=$(mktemp)
        while read -r transcript_id; do
            echo "transcript_id \"$transcript_id\""
        done < "$id_file" > "$temp_pattern_file"
        
        # Subset the GTF file using grep with the pattern file
        grep -F -f "$temp_pattern_file" "$merged_gtf" > "$output_gtf"
        
        # Clean up temporary file
        rm "$temp_pattern_file"
        
        # Count output lines
        output_count=$(wc -l < "$output_gtf")
        log_message "  -> Created $output_gtf with $output_count lines"
        log_message ""
    fi
done

log_message "All GTF subsetting completed successfully!"
log_message "Completed at: $(date)"
log_message "Log saved to: $LOG_FILE"