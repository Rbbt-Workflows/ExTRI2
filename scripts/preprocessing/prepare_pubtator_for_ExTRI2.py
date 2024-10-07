'''
This script processes multiple PubTator files, filters them based on a list of TF IDs,
and saves the filtered files to an output folder to use as input for the ExTRI2 pipeline.

Usage:
python prepare_pubtator_for_ExTRI2.py
'''

import os

def load_TF_ids(TF_list_file: str) -> set:
    '''
    Return set of TF ids from the TF list file.
    '''
    with open(TF_list_file, 'r') as f:
        return set(line.strip() for line in f)

def filter_pubtator_file(input_file: str, output_file: str, TF_ids: set) -> None:
    '''
    Filters a PubTator file to retain only the PMIDs with at least one gene annotation matching the valid genes.
    '''
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Separate the PubTator into a chunk for each PMID
        chunks = infile.read().strip().split('\n\n')

        for chunk in chunks:
            # Save into outfile if there's a TF ID in the annotations
            keep_pmid = False
            title, abstract, *entities = chunk.split('\n')

            for entity_line in entities:
                # Only check normalized entities
                entity_sections = entity_line.split('\t')
                if len(entity_sections) == 6:
                    pmid, start, end, mention, entity_type, normalization = entity_sections
                    if entity_type == 'Gene' and normalization in TF_ids:
                        keep_pmid = True

                elif len(entity_sections) > 6:
                    raise ValueError(f'Entity line has too many fields:\n{entity_line}')
        
            if keep_pmid:
                outfile.write(chunk + '\n\n')

def process_pubtator_files(input_folder: str, output_folder: str, TF_ids: set) -> None:
    '''
    Processes multiple PubTator files in the specified input folder, filters them, 
    and saves the results in the output folder.
    '''
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Iterate over each BioCXML.n.pubtator folder
    for n in range(10):
        subfolder = os.path.join(input_folder, f'BioCXML.{n}.pubtator/output')
        files = os.listdir(subfolder)
        total_files = len(files)

        for i, file in enumerate(files):
            input_file = os.path.join(subfolder, file)
            output_file = os.path.join(output_folder, f'{n}_{file}')

            # Discard PMIDs without TF IDs & save into output_file
            filter_pubtator_file(input_file, output_file, TF_ids)

            # Progress bar
            print(f"BioCXML.{n}: Processed {i+1}/{total_files} files", end='\r')
        print()

if __name__ == '__main__':

    DATA_FOLDER = '../../data/'
    gene_list_file = os.path.join(DATA_FOLDER, 'tf_entrez_code.list')
    input_folder   = os.path.join(DATA_FOLDER, 'tmp/BioCXML.pubtator')
    output_folder  = os.path.join(DATA_FOLDER, 'pubtator/')

    TF_ids = load_TF_ids(gene_list_file)
    process_pubtator_files(input_folder, output_folder, TF_ids)
