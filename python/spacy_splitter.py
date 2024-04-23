
import json
import sys
import spacy
from spacy.language import Language


#!pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.3/en_core_sci_md-0.5.3.tar.gz
nlp = spacy.load("en_core_sci_md")

# Ensure sentences are splitted at \n
@Language.component("set_custom_boundaries")
def set_custom_boundaries(doc):
    for token in doc[:-1]:
        if token.text == "\n":
            doc[token.i + 1].is_sent_start = True
    return doc
nlp.add_pipe("set_custom_boundaries", before="parser")


def split_texts(texts, chunk=300):
    '''
    Given a list of texts, return a list of [[[sent1, start, end],...]]
    '''
    # Pipe in chunks
    docs = []
    for i in range(0, len(texts), chunk):
        # Disable spacy components. Leave only 'parser', that separates text into sentences.
        docs  += list(nlp.pipe(texts[i:i+chunk], disable=['tagger', 'ner', 'lemmatizer', 'textcat']))

    # Convert to a list of lists for easy transfer to Ruby
    docs_sentences = []
    for doc in docs:
        docs_sentences.append(create_sentences_list(doc))
    
    return docs_sentences

    doc = nlp("" + text)
    sentences = [i for i in doc.sents]
    return sentences


def split_text(text):
    '''
    Given a text, return a list of [[sentence, start, end],...]
    '''
    # Disable spacy components. Leave only 'parser', that separates text into sentences.
    doc = nlp("" + text, disable=['tagger', 'ner', 'lemmatizer', 'textcat'])
    #sentences = [[sent.text, sent.start_char, sent.end_char] for sent in doc.sents]
    sentences = create_sentences_list(doc)

    return sentences


def create_sentences_list(doc):
    '''
    Given a nlp doc, return a list of sents in the format [[sentence, start, end],...], where
    2 sentences are kept separated only if the 1st ends with "[[.!?] ]|\n" AND the 2nd starts with uppercase.
    Returns:
        list
    '''
    sents = [i for i in doc.sents]  # Convert generator to a list for easier handling
    merged_sents = []
    temp_sent = ''
    start_char = 0

    for i, sent in enumerate(sents):
        
        temp_sent += sent.text

        # Separate to a new sentence if: 
        if  (i+1 == len(sents) or                             # It is the last sentence, or
             sents[i + 1].text.startswith('\n') or            # Next sentence starts with newline, or 
             sent.text.endswith('\n') or                      # This sentence ends with newline, or
             (sent.text.endswith(('.', '!', '?')) and         # ( Sentence ends with punctuation and
              sents[i + 1].start_char == sent.end_char + 1)): #   There is a space between sentences )

            assert doc.text[start_char:start_char+len(temp_sent)] == temp_sent, f"Sentence is not correctly aligned to the text"

            # Add sentence to doc
            end_char = sent.end_char
            merged_sents.append([temp_sent, start_char, end_char])

            # Reinitialise variables
            temp_sent = ""
            start_char = sents[i+1].start_char if i+1 != len(sents) else None
            
        else:
            # Add space between sentences if applicable
            temp_sent += ' '*(sents[i+1].start_char - sent.end_char)
    
    return merged_sents

if __name__ == "__main__":
    # Serialize the list to a JSON string
    json_input = sys.argv[1]
    texts = json.loads(json_input)
    spacy_transfer_list = split_texts(texts)
    print(json.dumps(spacy_transfer_list))