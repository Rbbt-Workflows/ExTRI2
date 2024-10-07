#!/bin/bash
# Source: https://www.ncbi.nlm.nih.gov/research/pubtator/api.html
# Given a list of PMIDs, return their annotated abstracts in pubtator format from PubTator3.
# Pubtator (still) only accepts GET commands, with a limit of 100 PMIDs per request.


# File containing the PubMed IDs + annotations
INPUT_PMIDS=$1 # All problematic PMIDs
OUTPUT_FILE=$2 # Output file where all results will be appended

# Base URL for the query
BASE_URL='https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/pubtator?pmids='

# Temporary file to hold chunks of 100 IDs
TMP_FILE="tmp_ids.txt"

# Split the file into chunks of 100 lines each
split -l 100 $INPUT_PMIDS $TMP_FILE

# Iterate over each chunk
for chunk in $(ls $TMP_FILE*); do
    echo "Processing file $chunk"
    
    # Read the IDs from the chunk and concatenate them with commas   
    ids=$(paste -sd, $chunk)
    
    # Construct the query URL
    query="${BASE_URL}${ids}"
    
    # Use wget to execute the query and append output to the single output file
    wget -qO- "$query" >> "$OUTPUT_FILE"
    #wget -qO- "$query"
    # wget gets stuck, so I just donwloaded it manually for now
    echo "$query"
    
    # Wait 1 second before making the next request to comply with rate limiting
    sleep 1

done

# Cleanup temporary files
rm $TMP_FILE**