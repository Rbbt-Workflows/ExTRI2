'''
Scripts surrounding manual normalisations
'''
import pandas as pd
import numpy as np
import itertools

# NFKB / AP1 specific functions 
def get_NFKB_masks(ExTRI2_df: pd.DataFrame, T='TF'):
    '''
    Return a mask of those entities to change to NFKB.
    T: either TF or TG.
    '''

    # Get the set of symbols that could be renormalized to NFKB
    NFKB_symbols = {'NFKB1', 'NFKB2', 'RELA', 'RELB'}

    # Filter for rows matching NFKB/REL symbols
    m_orig = ExTRI2_df[f'{T} Symbol'].str.upper().isin(NFKB_symbols)

    # Apply NFKB regex  only to the subset where m_nfkb is True
    col_regex = np.where(m_orig, apply_regex(ExTRI2_df[f'{T}'], 'NFKB'), '')

    # RULES:
    # 1. All whose Regex is NFKB will be changed to NFKB.
    m1 = col_regex == 'NFKB'

    # 2. Those normalized to NFKB1 that do not contain a 1 or P50 or NFKBIA will be also normalized to NFKB
    m2_1 = ExTRI2_df[f'{T} Symbol'].str.upper() == 'NFKB1'
    m2_2 = ~(ExTRI2_df[f'{T}'].str.contains('1') | ExTRI2_df[f'{T}'].str.upper().str.contains('P ?50'))
    m2 = m2_1 & m2_2

    m_final = m1 | m2
    return m_orig, m_final

def get_AP1_masks(ExTRI2_df: pd.DataFrame, T='TF'):
    '''
    T: either TF or TG.
    Returns:
      - a mask of all analyzed rows.
      - a mask of those entities that should be normalized to AP1.
    '''
    # Get the set of symbols containing FOS or JUN
    symbol_contains = ('FOS', 'JUN')
    dimer_symbols = set(ExTRI2_df[ExTRI2_df[f'{T} Symbol'].str.contains('|'.join(symbol_contains), case=False)][f'{T} Symbol'])

    # Filter for rows matching the symbols
    m_orig = ExTRI2_df[f'{T} Symbol'].isin(dimer_symbols)

    # Apply AP1 regex only to the subset where m_nfkb is True
    col_regex = np.where(m_orig, apply_regex(ExTRI2_df[f'{T}'], 'AP1'), '')

    # RULES:
    # Change to AP1 if it is AP1 or AP1 proteins
    m_final = np.isin(col_regex, ['AP1', 'AP1 PROTEINS'])
    return m_orig, m_final

def apply_regex(df: pd.Series, dimer: str) -> pd.Series:
    df = df.str.upper()
    df = df.str.replace(r"[-./;=\&{}()\[\]]", "", regex=True)

    if dimer == 'NFKB':
        df = df.str.replace(r'(NUCLE(AR|US))?( )*(TRANSCRIPTION(AL)? |REGULATORY )?FACTOR(S)?', 'NF', regex=True)
        df = df.str.replace(r'KAP(P?)( )?A( )?', 'K', regex=True)\
            .str.replace(r'KAP(P?)A( )?', 'K', regex=True)\
            .str.replace(r'BETA( )?', 'B', regex=True)
        df = df.str.replace(r'NF( )?(OF)?( )?(NF)?( )*', 'NF', regex=True)
        df = df.replace({r'REL A': 'RELA', r'REL B': 'RELB'}, regex=True)

    elif dimer == 'AP1':
        df = df.str.replace(r"ACTIVAT(ED|ION|ING|OR)( )?(OF )?PROTEIN", "AP", regex=True)
        df = df.str.replace(r"AP( AP)?( )?1", "AP1", regex=True)

    return df

def save_AP1_NFKB_breakdown_table(ExTRI2_df: pd.DataFrame, m_orig: pd.Series, m_final: pd.Series, T: str, dimer: str, config: dict) -> None:
    '''
    Create a table with the breakdown of sentences whose entity has been renormalized to AP1/NFKB
    
    Input:
        m_orig: mask of all sentences whose entity could potentially be renormalized
        m_final: mask of all sentences that are actually renormalized
        T: TF or TG
        dimer: AP1 or NFKB
    '''

    # Get paths & variables
    table_path = config['AP1_NFKB_breakdown_p']
    table_cols = config['AP1_NFKB_breakdown_cols']

    # Print number of sentences affected
    if T == 'TF':
        print(f"{m_final.sum()} rows ({m_final.sum()/len(m_final):.2%}) will have its {T} renormalized to {dimer}")
    elif T == 'TG':
        print(f"{m_final.sum()} rows ({m_final.sum()/len(m_final):.2%}) will be dropped as the {T} corresponds to {dimer}")

    # Initialise a list to store breakdown info on the renormalizations
    results = []

    # Iterate over each affected symbol
    symbols = set(ExTRI2_df[m_orig][f'{T} Symbol'].str.upper())

    for symbol in sorted(symbols):
        # Get sentences whose T is normalized to that symbol
        m_symbol = ExTRI2_df[f'{T} Symbol'].str.upper() == symbol

        # Initialise dictionary where to save results
        symbol_results = {}
        symbol_results['symbol'] = symbol
        symbol_results['entity'] = T
        symbol_results['dimer']  = dimer
        
        # Save percentage & regex names of unmodified and modified rows
        for state in 'unmodified', 'modified':
            # Get mask for each state
            m = m_symbol & ~m_final if state=='unmodified' else m_symbol & m_final

            # Use Regex to simplify & analyse the names more easily
            unique_x    = ExTRI2_df[m][T].drop_duplicates()
            regex_x     = ', '.join(sorted(apply_regex(unique_x, dimer).drop_duplicates()))
            #regex_x = highlight_symbol_in_text(regex_x, dimer)

            # Add to results dictionary
            symbol_results[f'% {state}']     = m.sum() / len(ExTRI2_df)
            symbol_results[f'regex {state}'] = regex_x

        results.append(symbol_results)

    # Convert the results list to a DataFrame & reorder columns
    results_df = pd.DataFrame(results)
    results_df = results_df[table_cols]

    # Append the results into the breakdown table
    table = pd.read_csv(table_path, sep='\t')
    table = pd.concat([table.astype(results_df.dtypes), results_df])
    table.to_csv(table_path, sep='\t', index=False)

    return

def correct_from_mask(ExTRI2_df, m, dimer):
    'Given a mask of sentences to change, correct them'
    # Change to dimer all those in mask. Normalise in the same way human and mouse/rat ones.
    ExTRI2_df.loc[m, "TF Symbol"]            = dimer
    ExTRI2_df.loc[m, "TF Id"]                = f'Complex:{dimer}'
    add_renormalisation_tag(ExTRI2_df, m, dimer)
    return ExTRI2_df

def correct_dimer(ExTRI2_df: pd.DataFrame, get_masks, dimer: str, config: dict, discarded_sents: list):
    for T in 'TF', 'TG':
        # Get masks for before & after renormalization
        m_orig, m_final = get_masks(ExTRI2_df, T)

        # Save the breakdown in a table
        save_AP1_NFKB_breakdown_table(ExTRI2_df, m_orig, m_final, T, dimer, config)

        # Renormalize TF entities, discard TG entities 
        if T == 'TF':
            ExTRI2_df = correct_from_mask(ExTRI2_df, m_final, dimer)        
        elif T == 'TG':
            discarded_sents.append(ExTRI2_df[m_final].copy())
            discarded_sents[-1]['Error'] = dimer
            ExTRI2_df = ExTRI2_df[~m_final]

    return ExTRI2_df

def fix_NFKB_AP1(ExTRI2_df: pd.DataFrame, config: dict) -> pd.DataFrame:
    '''Fix AP1 & NFKB normalisations'''

    # Create an empty table where to save the AP1/NFKB renormalization breakdown
    pd.DataFrame({v: [] for v in config['AP1_NFKB_breakdown_cols']}).to_csv(config['AP1_NFKB_breakdown_p'], sep='\t', index=False)
    
    # Initialise values
    discarded_sents = [] 
    ExTRI2_df['renormalisation'] = ''

    # Correct NFKB/AP1 entities & save results in the previous table   
    ExTRI2_df = correct_dimer(ExTRI2_df, get_NFKB_masks, 'NFKB', config, discarded_sents)
    ExTRI2_df = correct_dimer(ExTRI2_df, get_AP1_masks, 'AP1', config, discarded_sents)
    
    print(f"Breakdown by NCBI Symbol saved in {config['AP1_NFKB_breakdown_p']}")

    # Join & save discarded sentences
    pd.concat(discarded_sents).to_csv(config['NFKB_AP1_discarded_sents_p'], sep='\t')
    
    return ExTRI2_df


# General functions
def add_renormalisation_tag(ExTRI2_df, m, tag):
    ExTRI2_df.loc[m, 'renormalisation'] = np.where(
    ExTRI2_df.loc[m, 'renormalisation'] == '',
    tag,
    ExTRI2_df.loc[m, 'renormalisation'] + f';{tag}')
    return 

def renormalize(ExTRI2_df: pd.DataFrame, renormalized_sents_path: str) -> None:
    '''
    Renormalize common errors (excluding NFKB & AP1).
    Save all renormalized sentences in a table
    '''
    def apply_mappings(taxid_mapping:dict, mask: pd.Series):
        for taxid, (tg_symbol, tg_id) in taxid_mapping.items():
            taxid_mask = mask & (ExTRI2_df['TG TaxID'] == taxid)
            ExTRI2_df.loc[taxid_mask, "TG Symbol"] = tg_symbol
            ExTRI2_df.loc[taxid_mask, "TG Id"] = tg_id        

    def p21_to_CDKN1A(ExTRI2_df):
        '''Any instance of p21 must be normalized to CDKN1A'''
        mTG_p21 = (ExTRI2_df['TG'] == 'p21') & (ExTRI2_df['TG Symbol'].str.upper() != 'CDKN1A')
        print(f"{mTG_p21.sum()}\t{mTG_p21.sum()/len(mTG_p21):.2%}\tp21 is normalized to CDKN1A" )

        # Create a mapping of TaxID to (TG Symbol, TG Id)
        taxid_mapping = {
            '9606': ('CDKN1A', '1026'),
            '10090': ('Cdkn1a', '12575'),
            '10116': ('Cdkn1a', '114851')
        }
        
        # Renormalize Symbol & ID
        apply_mappings(taxid_mapping, mTG_p21)

        # Add tag
        add_renormalisation_tag(ExTRI2_df, mTG_p21, 'p21_to_CDKN1A')
        
        return

    def p53ps_to_TP53(ExTRI2_df):
        '''Correct/substitute all normalizations to “Trp53-ps” or “p53-ps” to the official gene symbol (Trp53 for mouse/rat, TP53 for human)'''
        mTG_p53 = (ExTRI2_df['TG'].str.upper() == 'P53') & (ExTRI2_df['TG Symbol'].str.contains('ps'))
        print(f"{mTG_p53.sum()}\t{mTG_p53.sum()/len(mTG_p53):.2%}\tp53-ps is normalized to its respective p53 symbol" )
        
        # Create a mapping of TaxID to (TG Symbol, TG Id)
        taxid_mapping = {
            '9606': ('TP53', '7157'),
            '10090': ('Trp53', '22059'),
            '10116': ('Tp53', '24842')
        }
        # Renormalize Symbol & ID
        apply_mappings(taxid_mapping, mTG_p53)

        # Add tag
        add_renormalisation_tag(ExTRI2_df, mTG_p53, 'p53-ps_to_p53')
        
        return

    print("Number of renormalized sentences and normalization:")
    p21_to_CDKN1A(ExTRI2_df)
    p53ps_to_TP53(ExTRI2_df)
    print()

    # Join & save discarded sentences
    if renormalized_sents_path is not None:
        ExTRI2_df[ExTRI2_df['renormalisation'] != ''].to_csv(renormalized_sents_path, sep='\t')

    return

def discard(ExTRI2_df: pd.DataFrame, discarded_sents_path: str) -> pd.DataFrame:
    '''
    Discard sentences with common errors (excluding NFKB & AP1).
    Save all discarded sentences in a table
    '''

    discarded_sents = []
    num_sents = len(ExTRI2_df)

    def discard_NAT(ExTRI2_df):
        '''long noncoding antisense RNAs (NATs) mislabelled as TFs'''
        # They are recognizeable by the suffix -AS[1-3]
        m1 = ExTRI2_df['TF'].str.contains('-AS[1-3]')
        m2 = ExTRI2_df['Sentence'].str.contains('\[TF\]-AS[1-3]')
        m_NAT = m1 | m2
        print(f"{m_NAT.sum()}\t{m_NAT.sum()/num_sents:.2%}\tTheir TF contains -AS[1-3]" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m_NAT].copy())
        discarded_sents[-1]['Error'] = 'lnc antisense RNAs'
        
        return ExTRI2_df[~m_NAT]

    def discard_circRNAs(ExTRI2_df):
        '''circRNAs misnormalized as TFs must be discarded'''
        m1 = ExTRI2_df['TF'].str.startswith('circ')
        m2 = ExTRI2_df['Sentence'].str.contains('circ(?:RNAs?)?[ -]?\[TF\]')
        m_circ = m1 | m2
        print(f"{m_circ.sum()}\t{m_circ.sum()/num_sents:.2%}\tTheir TF are circRNAs" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m_circ].copy())
        discarded_sents[-1]['Error'] = 'circRNA'
        
        return ExTRI2_df[~m_circ]

    def discard_NLRP3_inflammasome(ExTRI2_df):
        '''NLRP3 inflammasome normalised to NLRP3 must always be discarded'''
        for T in ('TF', 'TG'):
            m_NRLP3 = (ExTRI2_df[f'{T}'].str.upper() == 'NLRP3') & (ExTRI2_df['Sentence'].str.contains(f'\[{T}\] inflammasome'))
            print(f"{m_NRLP3.sum()}\t{m_NRLP3.sum()/num_sents:.2%}\tTheir {T} (NLRP3) is followed by inflammasome but normalised to NLRP3" )

            # Add discarded sentences to list
            discarded_sents.append(ExTRI2_df[m_NRLP3].copy())
            discarded_sents[-1]['Error'] = 'NLRP3 inflammasome'

            ExTRI2_df = ExTRI2_df[~m_NRLP3]
        
        return ExTRI2_df

    def discard_fusion_genes(ExTRI2_df):
        '''Discard Common fusion genes misnormalized as TGs'''

        # TODO - This is incomplete. There's way more fusion genes that should be discarded apart from these ones.

        fusion_genes = [('ABL1', 'BCR'), ('FLI1','EWSR1')]
        pairs = [';'.join(p) for pair in fusion_genes for p in itertools.permutations(pair, 2) ]
        m_fusion = ExTRI2_df['TG Symbol'].isin(pairs)
        print(f"{m_fusion.sum()}\t{m_fusion.sum()/num_sents:.2%}\tThe TG is a fusion gene" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m_fusion].copy())
        discarded_sents[-1]['Error'] = 'fusion gene'

        return ExTRI2_df[~m_fusion]

    def discard_invalid_TF_function(ExTRI2_df):
        '''Discard sentences with mentions invalid of TG function'''
        m_funct = ExTRI2_df['Sentence'].str.contains(r'\[TG] (?:pathway|signall?ing|axis|program)')
        print(f"{m_funct.sum()}\t{m_funct.sum()/num_sents:.2%}\tTheir TG is followed by pathway/signalling/axis/program" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m_funct].copy())
        discarded_sents[-1]['Error'] = 'invalid function'

        return ExTRI2_df[~m_funct]

    def discard_MDM2_TP53(ExTRI2_df):
        '''Sentences with MDM2-TP53 pairs must be removed: they're always a PPI'''
        m = ExTRI2_df['TF Symbol'].str.upper().str.contains('MDM2')
        m &= ExTRI2_df['TG Symbol'].str.upper().str.contains('P53')
        print(f"{m.sum()}\t{m.sum()/num_sents:.2%}\tMDM2-TP53 pair, which is always a PPI" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m].copy())
        discarded_sents[-1]['Error'] = 'MDM2-TP53'

        return ExTRI2_df[~m]

    def discard_CD_cells(ExTRI2_df):
        '''If TF/TG is CD4, CD8A, CD8B, CD74, CD34, and followed by +/positive, discard'''
        # This only returns 14 cases, indicating that the model is already mostly good at
        # identifying that CDX+ refers to a cell type and is therefore an incorrect sentence.

        for T in ['TF', 'TG']:

            m_base = ExTRI2_df[T].str.contains(r'^CD(?:4|8A|8B|74|34)(?!\d)')
            m1 = ExTRI2_df['Sentence'].str.contains(f'\[{T}\]\ ?[\[\()]?\+')
            m2 = ExTRI2_df['Sentence'].str.contains(f'\[{T}\] ?positive')
            m = m_base & (m1 | m2)
            print(f"{m.sum()}\t{m.sum()/num_sents:.2%}\t{T} is CD(4|8A|8B|74|34) positive (indicative of a cell, not a gene)" )

            # Add discarded sentences to list
            discarded_sents.append(ExTRI2_df[m].copy())
            discarded_sents[-1]['Error'] = f'{T} is CD cell'

            ExTRI2_df = ExTRI2_df[~m]
        
        return ExTRI2_df

    def discard_of_entity(ExTRI2_df):
        m = ExTRI2_df['TG'] == 'of'
        m |= ExTRI2_df['TF'] == 'of'

        print(f"{m.sum()}\t{m.sum()/num_sents:.2%}\tPreposition of incorrectly identified as a gene" )

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m].copy())
        discarded_sents[-1]['Error'] = 'of_preposition'

        return ExTRI2_df[~m]
    
    def discard_translation(ExTRI2_df):
        '''Most sentences containing translation are wrong, so they'll be tagged'''
        m = ExTRI2_df['Sentence'].str.lower().str.contains('translat')
        print(f"{m.sum()}\t{m.sum()/num_sents:.2%}\tThey contain translation in them (found to only be correct 40% of the time)")

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m].copy())
        discarded_sents[-1]['Error'] = f'Translation'

        return ExTRI2_df[~m]

    def discard_autoregulation(ExTRI2_df):
        '''Most sentences showing autoregulation are wrong, so they'll be tagged'''
        m = ExTRI2_df['TF Symbol'].str.upper() == ExTRI2_df['TG Symbol'].str.upper()
        print(f"{m.sum()}\t{m.sum()/num_sents:.2%}\tThey are autoregulations (found to only be correct ~10% of the time)")

        # Add discarded sentences to list
        discarded_sents.append(ExTRI2_df[m].copy())
        discarded_sents[-1]['Error'] = f'Autoregulation'

        return ExTRI2_df[~m]

    print(f"Number of discarded sentences and percentage from total ({num_sents} sentences) and reasoning:")
    ExTRI2_df = discard_NAT(ExTRI2_df)
    ExTRI2_df = discard_circRNAs(ExTRI2_df)
    ExTRI2_df = discard_NLRP3_inflammasome(ExTRI2_df)
    ExTRI2_df = discard_fusion_genes(ExTRI2_df)
    ExTRI2_df = discard_invalid_TF_function(ExTRI2_df)
    ExTRI2_df = discard_MDM2_TP53(ExTRI2_df)
    ExTRI2_df = discard_CD_cells(ExTRI2_df)
    ExTRI2_df = discard_of_entity(ExTRI2_df)
    ExTRI2_df = discard_translation(ExTRI2_df)
    ExTRI2_df = discard_autoregulation(ExTRI2_df)
    
    print()

    # Join & save discarded sentences
    if discarded_sents_path is not None:
        pd.concat(discarded_sents).drop(columns=['TF_type']).to_csv(discarded_sents_path, sep='\t')

    return ExTRI2_df


