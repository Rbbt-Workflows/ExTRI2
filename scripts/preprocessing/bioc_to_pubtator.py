'''
Convert a truncated bioc file into a pubtator file.
Used inside bioc_to_pubtator.sh
'''

#!pip install bconv
import bconv
import sys
import os

def coll_to_pubtator_str(coll, filename: str, problematic_pmids_file: str) -> str:
    '''Convert the title and abstract, and its entities, of the collection into a PubTator format string'''

    pubtator_str = ''
    for doc in coll:

        # TEST
        # If title is not present, raise error. It should be present in all documents.
        assert (doc[0].type == 'front' or doc[0].type == 'title'), f"{doc.id} doesn't start with a title"
        # DEBUG: print(doc.id, [doc[i].type for i in range(len(doc))], file=sys.stderr)
        
        # If a section is not abstract, ignore PMID.
        if (len(doc)>1) and (any([('abstract' not in section.type) for section in doc[1:]])):
            with open(problematic_pmids_file, 'a') as f:
                    print(f"{filename}:{doc.id}:({len(doc)}){set([section.type for section in doc[1:]])}", file=sys.stderr)
                    f.write(f'{filename}:{doc.id}:NotAbstractSections:({len(doc)}){[section.type for section in doc[1:]]}\n') 
            continue
        

        # WRITE TITLE
        # Add spaces at the end if title.end != abstract.start+1
        pubtator_str += doc.id + '|t|' + doc[0].text
        pubtator_str += ' '*(doc[1].start - doc[0].end - 1)+'\n' if len(doc) > 1 else '\n'


        # WRITE ABSTRACT
        # Write [first section of the] abstract
        pubtator_str += doc.id + '|a|'
        pubtator_str += doc[1].text if len(doc) > 1 else ''

        # If more than 1 abstract section, separate with |
        if len(doc)>2:
            for i, section in enumerate(doc[2:], start=2):
                # Add section preceded by spaces if prev_section.end != section.start+1
                pubtator_str += ' '*(section.start - doc[i-1].end - 1) + '|' + section.text

        pubtator_str += '\n'

        # ADD ANNOTATIONS following PubTator format
        gene_count = 0 # Count gene mentions
        for i, section in enumerate(doc):          
            for sentence in section:
                for ent in sentence.entities:
                    metadata = f'{ent.metadata["type"]}'
                    metadata += f'\t{ent.metadata["identifier"]}' if 'identifier' in ent.metadata else ''
                    pubtator_str += f'{doc.id}\t{ent.start}\t{ent.end}\t{ent.text}\t{metadata}\n'

                    # Count the number of gene mentions
                    if ent.metadata['type'] == 'Gene':
                        gene_count += 1


        # Add blank line
        pubtator_str += '\n'

        # Discard if no gene mention
        if gene_count < 1:
            pubtator_str = ''

    return pubtator_str

def main(filename: str, input_bioc_file: str, output_pubtator_file: str, problematic_pmids_file: str) -> None:
    '''Convert the BioC into PubTator, or return the problematic PMID'''
    pubtator_str = ''

    # Convert the truncated BioC to a pubtator
    try:
        coll = bconv.load(input_bioc_file, fmt='bioc_xml')
    except AssertionError as e:
        print(f"{filename}:AssertionError:{str(e).split(':')[0]}", file=sys.stderr)
        
        # Extract the PMID that causes the error
        error_PMID = str(e).split(':')[0].split()[-1]
        assert error_PMID.isdigit(), f"ERROR: PMID extracted from AssertionError ({error_PMID}) is not an integer"

        # Print error_PMID in the standard output to be removed from the collection
        print(error_PMID)
        sys.exit(1) # Exit with error code to indicate an issue

    pubtator_str = coll_to_pubtator_str(coll, filename, problematic_pmids_file)
    with open(output_pubtator_file, 'a') as f:
        f.write(pubtator_str)


if __name__ == '__main__':
    # filename, input_bioc_file, output_pubtator_file, problematic_pmids_file
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
