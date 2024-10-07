#!/bin/bash
# Convert BioC to Pubtator, only keeping title and abstract from those PMIDs with Gene mentions
# If abstract has more than 1 section (e.g. abstract_title1, abstract, etc.), separate by |
# If PMID cannot be processed, obtain from PubTator through POST
#
# Usage: bioc_to_pubtator.sh <input_dir> <output_dir>

# Define truncate function
truncate_bioc() {
    
    local INPUT_FILE=$1
    local OUTPUT_FILE=$2

    # In PMIDs with abstract, get title & abstract
    grep -oP '^.*<infon key="type">abstract.*?<\/passage>' $INPUT_FILE | sed 's/$/<\/document>/' > $OUTPUT_FILE
    # In PMIDs without abstract, get only title
    grep -v '<infon key="type">abstract' $INPUT_FILE | grep -oP '^.*?<infon key="section_type">TITLE<\/infon>.*?</passage>' | sed 's/$/<\/document>/' >> $OUTPUT_FILE

    # Only keep docs that have Gene annotations
    grep 'key="type">Gene</infon>' "$OUTPUT_FILE" > temp_file && mv temp_file "$OUTPUT_FILE"
    # Add head and tail
    {
        head -1 "$INPUT_FILE"
        cat "$OUTPUT_FILE"
        tail -1 "$INPUT_FILE"
    } > temp_file && mv temp_file "$OUTPUT_FILE"

    # Save ignored docs IDs (no abstract nor title) in the IGNORED_DOCS file. 
    grep -v '<infon key="type">abstract' $INPUT_FILE | grep -v '<infon key="section_type">TITLE<\/infon>' | grep -Pv '<\/?collection>' | grep -o '<id>.*</id>' >> $IGNORED_DOCS_FILE

}


# Activate python environment
source .bioc/bin/activate

# Get input and output directories from stdin
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi
INPUT_BIOC_DIR=$1
OUTPUT_DIR=$2


# Define variables
OUTPUT_PUBTATOR_DIR=$OUTPUT_DIR/output
PROBLEMATIC_PMIDS_FILE=$OUTPUT_DIR/problematic_PMIDs
IGNORED_DOCS_FILE=$OUTPUT_DIR/ignored_docs
NUM_FILES=$(ls "$INPUT_BIOC_DIR" | wc -l)
CURRENT_FILE=0
OUTPUT_PUBTATOR_FILE=$OUTPUT_PUBTATOR_DIR/$CURRENT_FILE.pubtator

# Initialise directories and variables
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_PUBTATOR_DIR
>$IGNORED_DOCS_FILE
>$PROBLEMATIC_PMIDS_FILE


# Process each BioC collection in the input directory
for INPUT_FILE in "$INPUT_BIOC_DIR"/*; do

    FILE_NAME=$(basename "$INPUT_FILE")

    # Truncate file to contain only title and abstract & keep only PMIDs with gene mentions
    truncate_bioc $INPUT_FILE truncated_bioc

    # Convert truncated bioc into pubtator, ignoring problematic pmids
    while : ; do

        # Convert truncated bioc into pubtator. Raise error and save problematic PMID if AssertionError
        error_PMID=$(python bioc_to_pubtator.py $FILE_NAME truncated_bioc $OUTPUT_PUBTATOR_FILE $PROBLEMATIC_PMIDS_FILE)
        exit_status=$?

        if [ $exit_status -eq 0 ]; then
            # Conversion succeeded
            break
        elif [ $exit_status -eq 1 ]; then
            # If AssertionError, remove problematic PMID from collection and save it into PROBLEMATIC_PMIDS_FILE
            echo "$FILE_NAME: Processing failed due to PMID: $error_PMID. Removing and trying again."
            grep -v "$error_PMID" truncated_bioc > temp_file && mv temp_file truncated_bioc
            echo "$FILE_NAME:$error_PMID:AssertionError" >> $PROBLEMATIC_PMIDS_FILE
        else
            echo "An unknown error occurred: " $exit_status
        fi

    done
    
    # Every 100 collections, start a new pubtator file
    if ((CURRENT_FILE % 100 == 0)) || ((CURRENT_FILE + 1 == NUM_FILES)); then
        OUTPUT_PUBTATOR_FILE=$OUTPUT_PUBTATOR_DIR/$((CURRENT_FILE+100)).pubtator

        # Progress bar update
        PERCENT=$((CURRENT_FILE * 100 / NUM_FILES))
        BAR=$(printf '%*s' $((PERCENT/2)) | tr ' ' '#')
        printf "\rProgress: [%-50s] %d%% (%d/%d)" "$BAR" "$PERCENT" "$CURRENT_FILE" "$NUM_FILES"

    fi
    CURRENT_FILE=$((CURRENT_FILE+1))

done

rm truncated_bioc

echo
echo Number of problematic docs:  $(cat $IGNORED_DOCS_FILE | wc -l)

echo Number of problematic PMIDs: $(cat $PROBLEMATIC_PMIDS_FILE | wc -l)
echo Obtaining them from PubTator through GET...

# Get problematic PMIDs from the file. Remove PMC ones (they are not in PubTator format)
awk '{split($1, arr, ":"); print arr[2]}' $PROBLEMATIC_PMIDS_FILE | grep -v "PMC" > tmp_PMIDs
bash get_pmid_list.sh tmp_PMIDs $OUTPUT_PUBTATOR_FILE
echo Done
