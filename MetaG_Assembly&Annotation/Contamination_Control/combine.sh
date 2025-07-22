#!/bin/bash

# Usage: ./combine_fasta_with_unique_headers.sh /path/to/fasta_dir output.fa

set -e

INPUT_DIR="$1"
OUTPUT_FILE="$2"

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_FILE" ]]; then
    echo "Usage: $0 /path/to/fasta_files/ combined_output.fa"
    exit 1
fi

# Empty the output file if it exists
> "$OUTPUT_FILE"

# Find all fasta files (both compressed and uncompressed)
find "$INPUT_DIR" -type f \( -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" -o -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fna.gz" \) | while read -r fasta; do

    base=$(basename "$fasta")
    prefix="${base%%.*}"  # Extract everything up to the first '.'

    echo "Processing $base with prefix $prefix..."

    if [[ "$fasta" == *.gz ]]; then
        zcat "$fasta"
    else
        cat "$fasta"
    fi | awk -v prefix="$prefix" '
        BEGIN { OFS = "" }
        /^>/ {
            header = substr($0, 2)
            print ">" prefix "_" header
            next
        }
        { print }
    ' >> "$OUTPUT_FILE"

done

echo "✅ Combined FASTA written to: $OUTPUT_FILE"
