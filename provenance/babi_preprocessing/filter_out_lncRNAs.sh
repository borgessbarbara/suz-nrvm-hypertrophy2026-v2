#!/bin/bash

DIR="/home/vmaracaja/shared-projects/sebauz_sr"
cd "$DIR"

idents=("CNE_24_HISAT2" "CNE_24_STAR" "CNE_6_24_HISAT2" "CNE_6_24_STAR")

for ident in "${idents[@]}"; do
    tx_file="${ident}.annotated.lncRNA_tx.gtf"
    id_list="${ident}_transcript_ids.txt"
    merged_file="merged_${ident}.annotated.gtf"
    output_file="filtered_merged_${ident}.annotated.gtf"

    # Extract transcript_ids
    awk -F 'transcript_id "' '{if (NF>1) {split($2,a,"\""); print a[1]}}' "$tx_file" | sort | uniq > "$id_list"

    # Remove lines from merged file with those transcript_ids
    awk '
        BEGIN {
            while ((getline id < "'$id_list'") > 0) {
                ids[id] = 1
            }
            close("'$id_list'")
        }
        {
            match($0, /transcript_id "([^"]+)"/, arr)
            if (arr[1] == "" || !(arr[1] in ids)) print $0
        }
    ' "$merged_file" > "$output_file"

    echo "Filtered $merged_file -> $output_file"
done