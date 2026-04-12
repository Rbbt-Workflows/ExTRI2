#!/bin/bash
# Create a file with all Gene mentions from PubTator3

INPUT_FOLDER=data/preprocessing/BioCXML.pubtator

OUTPUT_FILE=data/preprocessing/pubtator_gene_mentions.txt
GENE_IDS=data/preprocessing/pubtator_gene_ids.txt

TMP_FILE=data/preprocessing/tmp

> $OUTPUT_FILE
> $GENE_IDS
> $TMP_FILE


# Iterate over all subfolders of input_folder
for subfolder in $INPUT_FOLDER/*; do
    echo "Processing $subfolder"
    for file in $subfolder/output/*; do
        # Extract gene mentions and append to tmp file
        grep -E "\sGene\s*[0-9]" "$file" | awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $4, $5, $6}' | sort | uniq  >> "$OUTPUT_FILE"
    done

    # Only save Gene IDs
    awk -F'\t' '$1 <= 38776000 {print $4}' "$OUTPUT_FILE" | sort | uniq >> "$TMP_FILE"
done

# Get unique gene IDs
cat $TMP_FILE | sort | uniq > $GENE_IDS
rm $TMP_FILE