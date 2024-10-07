import pandas as pd
import math

TaxID_to_species = {
    '9606': 'human',
    '10090': 'mouse',
    '10116': 'rat'
}

# CONFIG
def load_config() -> dict:
    config = {}

    # Folders
    RESULTS_P       = '../../results/'
    DATA_P          = '../../data/'
    CLASSIFIERS_P   = '../../classifiers_training/data/'

    TEMP_P          = DATA_P + 'tmp/'
    VALIDATION_P    = DATA_P + 'validation/'
    TO_VALIDATE_P   = VALIDATION_P + 'to_validate/'
    VALIDATED_P     = VALIDATION_P + 'validated/'
    config['TO_VALIDATE_P'] = TO_VALIDATE_P



    # Results
    config['raw_ExTRI2.gz_p']          = RESULTS_P + 'ExTRI2.tsv.gz'
    config['final_ExTRI2_p']        = RESULTS_P + 'ExTRI2_final_resource.tsv'
    config['final_validated_p']     = RESULTS_P + 'validated_sentences.tsv'

    # Temporary files
    config['nonvalid_sample_p']     = TEMP_P + 'ExTRI2_nonvalid_subset.tsv'
    config['valid_pre_renorm_p']    = TEMP_P + 'ExTRI2_valid_prenorm.tsv'         # Used for validation purposes


    config['tf_ids_p']              = DATA_P + 'tf_entrez_code.list'
    config['training_data_p']       = CLASSIFIERS_P + 'TRI_data.tsv'

    # Validated
    config['stage1_validated_p']             = VALIDATED_P + '550_Iter1_validated_sentences_recheck.txt'
    config['stage2_validated_p']             = VALIDATED_P + '07_24_sents_to_validate_AL.txt'
    config['stage3_validated_p']             = VALIDATED_P + '07_28_sents_to_validate_3_AL_v2.txt'
    config['stage4_coTF_validated_p']        = VALIDATED_P + '19_08_coTF_to_validate_2_AL.txt'
    config['stage4_balanced_dbTF_validated_p'] = VALIDATED_P + '19_08_balanced_dbTF_to_validate_AL_v5_complete.txt'
    config['stage5_coTF_validated_p']        = VALIDATED_P + '26_08_likely_coTF_to_validate_AL.txt'
    config['stage6_coTF_validated_p']        = VALIDATED_P + '28_08_likely_coTF_to_validate_AL.txt'
    config['stage6_ll_coTF_validated_p']     = VALIDATED_P + '28_08_unlikely_coTF_to_validate_AL.txt'

    # TODO - Modify this to the real path
    config['stage7_dbTF_validated_p']        = TO_VALIDATE_P + '01_10_dbTF_to_validate.tsv'
    config['stage7_coTF_validated_p']        = TO_VALIDATE_P + '01_10_coTF_to_validate.tsv'
    config['stage7_ll_coTF_validated_p']     = TO_VALIDATE_P + '01_10_ll_coTF_to_validate.tsv'


    # To validate
    # TODO - Add how were those first 3 created.
    config['stage1_dbTF_to_validate_p']          = TO_VALIDATE_P + '04_28_sentences_to_validate.tsv'
    config['stage2_dbTF_to_validate_p']          = TO_VALIDATE_P + '07_24_sentences_to_validate_2.tsv'
    config['stage3_dbTF_to_validate_p']          = TO_VALIDATE_P + '07_28_sentences_to_validate_3.tsv'
    config['stage3_GTF_to_validate_p']           = TO_VALIDATE_P + '07_28_sentences_to_validate_3_GTFs.tsv'

    config['stage4_balanced_dbTF_to_validate_p'] = TO_VALIDATE_P + '19_08_balanced_dbTF_to_validate.tsv'
    config['stage4_random_coTF_to_validate_p']   = TO_VALIDATE_P + '19_08_random_coTF_to_validate.tsv'
    config['stage4_balanced_coTF_to_validate_p'] = TO_VALIDATE_P + '19_08_balanced_coTF_to_validate.tsv'

    config['stage5_likely_coTF_to_validate_p']   = TO_VALIDATE_P + '26_08_likely_coTF_to_validate.tsv'
    config['stage5_unlikely_coTF_to_validate_p'] = TO_VALIDATE_P + '26_08_unlikely_coTF_to_validate.tsv'

    config['stage6_likely_coTF_to_validate_p']   = TO_VALIDATE_P + '28_08_likely_coTF_to_validate.tsv'
    config['stage6_unlikely_coTF_to_validate_p'] = TO_VALIDATE_P + '28_08_unlikely_coTF_to_validate.tsv'

    config['stage7_dbTF_to_validate_p']          = TO_VALIDATE_P + '01_10_dbTF_to_validate.tsv'
    config['stage7_coTF_to_validate_p']          = TO_VALIDATE_P + '01_10_coTF_to_validate.tsv'
    config['stage7_ll_coTF_to_validate_p']       = TO_VALIDATE_P + '01_10_ll_coTF_to_validate.tsv'

    return config


# GET FINAL JOINED DATASET
def join_all_validated_dfs(config) -> pd.DataFrame:
    '''
    Return a joined df of all validated dfs, as there are some differences between their column names & values.
    '''

    def load_tsv(path: str) -> pd.DataFrame:
        return pd.read_csv(config[path], sep='\t', header=0, dtype=str)

    validated_dfs = {
        's1':         load_tsv('stage1_validated_p'),
        's2':         load_tsv('stage2_validated_p'),
        's3':         load_tsv('stage3_validated_p'),
        's4_dbTF':    load_tsv('stage4_balanced_dbTF_validated_p'),
        's4_coTF':    load_tsv('stage4_coTF_validated_p'),
        's5_coTF':    load_tsv('stage5_coTF_validated_p'),
        's6_coTF':    load_tsv('stage6_coTF_validated_p'),
        's6_coTF_ll': load_tsv('stage6_ll_coTF_validated_p'),
        's7_dbTF':    load_tsv('stage7_dbTF_validated_p'),
        's7_coTF':    load_tsv('stage7_coTF_validated_p'),
        's7_coTF_ll': load_tsv('stage7_ll_coTF_validated_p'),
    }

    # Drop empty rows
    for key in validated_dfs.keys():
        validated_dfs[key].dropna(how='all', inplace=True)

    # Drop non-validated rows
    dropped_rows = []
    for stage in ('s4_coTF', 's5_coTF'):
        m = validated_dfs[stage]['Valid?'].isna()
        validated_dfs[stage] = validated_dfs[stage][~m]

        dropped_rows.append(m.sum())
    print(f"{dropped_rows[0]} and {dropped_rows[1]} rows have been removed from the stage 4 & 5 validated datasets respectively as they were not validated.")
    # This drop was not done when creating the random & balanced datasets of stage 6. This has caused us to have to validate additional balanced dbs (stage 7).

    # Add missing columns: TF type
    for stage in ('s1', 's2', 's4_dbTF'):
        validated_dfs[stage]['TF_type'] = 'dbTF'
    for stage in ('s4_coTF', 's5_coTF', 's6_coTF'):
        validated_dfs[stage]['TF_type'] = 'coTF'
    validated_dfs['s6_coTF_ll']['TF_type'] = 'll_coTF'

    # Add missing columns: method
    for stage in ('s1', 's2', 's3'):
        validated_dfs[stage]['method'] = 'random'
    validated_dfs['s4_dbTF']['method'] = 'balanced'

    # Add missing columns: Sentence
    validated_dfs['s1']['Sentence'] = validated_dfs['s1']['span_sentence'].str.replace(r"<TF>.*<\/TF>", "[TF]", regex=True).str.replace(r"<TG>.*<\/TG>", "[TG]", regex=True)

    # Add stage column (except for stage 7, for which stage is defined below)
    validated_dfs = {stage: (df if 's7' in stage else df.assign(stage=stage)) for stage, df in validated_dfs.items()}

    # Rename columns
    for stage in ('s7_dbTF', 's7_coTF', 's7_coTF_ll'):
        validated_dfs[stage].drop(columns=['TF_type'], inplace=True)
        validated_dfs[stage].rename(columns={'val TF_type': 'TF_type'}, inplace=True)     
        validated_dfs[stage]['stage'] = validated_dfs[stage]['stage'] + '_' + validated_dfs[stage]['TF_type']

        # TODO - ERASE
        validated_dfs[stage]['Valid?'] = 'T'

    col_renamings = {'dbTF_norm': 'TF Symbol', 'Gene norm': 'TG Symbol', 'Gene': 'TG',
                 'Incorrect?': 'TF_is_incorrect', 'Incorrect?.1': 'TG_is_incorrect',
                 'Mention': 'TF_correct_mention', 'Mention.1': 'TG_correct_mention',
                 "Incorrect? TF": "TF_is_incorrect", "Incorrect? TG": "TG_is_incorrect",
                 "Correct name TF": "TF_correct_mention", "Correct name TG": "TG_correct_mention",
                 'Correct name': 'TF_correct_mention', 'Correct name.1': 'TG_correct_mention',
                 'Species': 'TF Species', 'TF type': 'TF_type',
                 }
    
    for stage, df in validated_dfs.items():
        df.rename(columns={old: new for old, new in col_renamings.items() if old in df.columns}, inplace=True)

    # Join in one dataframe
    validated_df = pd.concat(validated_dfs.values(), axis=0)

    # Drop sentences with unexpected missing values that shouldn't contain missing values
    drop_mask = validated_df['TG Id'].isna() | validated_df['Valid?'].isna() | ((validated_df['Valid?'] == 'F') & validated_df['true_label'].isna())
    validated_df = validated_df[~drop_mask]
    print(f"{drop_mask.sum()} sentences have been dropped because they contained missing values.")

    # Make all #SentenceID have the same form: PMID:Sent:TF:TG & define their PMID
    m_v = validated_df['#SentenceID'].str.startswith('PMID')
    validated_df.loc[m_v, '#SentenceID'] = validated_df[m_v]['#SentenceID'].apply(lambda row: ":".join(row.split(':')[i] for i in [1,4,5,6]))
    validated_df['PMID'] = validated_df['#SentenceID'].apply(lambda row: row.split(':')[0])

    # Rename MoR & Label
    validated_df['MoR'] = validated_df['MoR'].replace('REPRESION', 'REPRESSION')
    validated_df['true_MoR'] = validated_df['true_MoR'].replace({'repRESSION': 'REPRESSION', 'repression': 'REPRESSION'})

    # Join other_issues column
    cols_to_join = ['OI2', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 21', 'Unnamed: 25']
    validated_df['Other issues'] = validated_df[['Other issues'] + cols_to_join].apply(lambda x: ';'.join(sorted(x.dropna())), axis=1)
    validated_df.drop(columns=cols_to_join, inplace=True)

    # Assertions
    assert validated_df['method'].isna().sum() == 0, 'Method is not defined in some stages'
    assert (~validated_df['TF_type'].isin(['dbTF', 'GTF', 'coTF', 'll_coTF'])).sum() == 0, 'TF type category error. It should be dbTF, coTF, or ll_coTF.'

    return validated_df

def join_validated_df_with_valid_df(validated_df, valid_df, state='prerenorm'):
    '''
    Merge the validated_df with the valid_df. Used both for
    Return the merged and not merged sentences
    '''

    # Do some changes on valid_df so it has the same format as validated_df
    valid_df = valid_df.rename(columns={'Valid': 'Label'})
    valid_df['Label'] = valid_df['Label'].replace({'Valid': 'TRUE', 'Non valid': 'FALSE'})

    # Only keep the relevant columns from validated_df
    validated_df_relevant_cols = ['#SentenceID', 'Sentence', 'TF_type',
        'Label', 'MoR', 'Valid?', 'dir-*', 'true_label', 'true_MoR', 
        'Other issues', 'Explanation', 'method', 'stage',
        'TF Id', 'TF Symbol', 'TF_is_incorrect', 'TF_correct_mention',
        'TG Id', 'TG Symbol', 'TG_is_incorrect', 'TG_correct_mention', 
    ]

    validated_df = validated_df[validated_df_relevant_cols].copy()
    
    # Merge the validated_df with the valid_df_{state} df
    merged_df = pd.merge(validated_df, valid_df, on=['#SentenceID', 'Sentence'], how='inner', suffixes=('_validated', f'_{state}'))
    print(f"{len(merged_df)}/{len(validated_df)} rows are merged.")

    # Get non_merged columns
    non_merged_df = validated_df[~validated_df['#SentenceID'].isin(merged_df['#SentenceID'])]

    return merged_df, non_merged_df

def find_remaining_validated_in_candidate_sents(non_merged_df, config, chunk_size=1_000_000):
    '''Merge the sentences not found in valid_df with ExTRI2_df more efficiently'''

    candidate_sents = pd.read_csv(config['raw_ExTRI2.gz_p'], sep='\t', header=1, chunksize=chunk_size, keep_default_na=False)
    
    # Create empty list to collect matches
    found_sents_list = []

    for i, chunk in enumerate(candidate_sents):
        # Filter the chunk by 'Sentence' to reduce the number of rows we process further
        chunk_filtered = chunk[chunk['Text'].isin(non_merged_df['Sentence'])].copy()

        # Postprocess only the filtered chunk
        if not chunk_filtered.empty:
            chunk_filtered['#SentenceID'] = chunk_filtered['#SentenceID'].apply(lambda row: ":".join(row.split(':')[i] for i in [1, 4, 5, 6]))

            # Find rows where both '#SentenceID' and 'Sentence' match
            m = (chunk_filtered['#SentenceID'].isin(non_merged_df['#SentenceID'])) & (chunk_filtered['Text'].isin(non_merged_df['Sentence']))
            if m.sum() > 0:
                found_sents_list.append(chunk_filtered[m])


        print(f'Processed {(i+1)*chunk_size} rows from the candidate sentences', end='\r')

    # Changes to make candidate_sents format + validated sents format
    found_sents = pd.concat(found_sents_list, axis=0)
    found_sents = found_sents.rename(columns={'Text': 'Sentence', 'Valid': 'Label', 'Gene': 'TG', 'IG Id': 'TG Id'})
    found_sents['Label'] = found_sents['Label'].replace({'Valid': 'TRUE', 'Non valid': 'FALSE'})

    # Change column names for consistency with the previous merged df (as candidate_sents don't have this metadata)
    non_merged_df = non_merged_df.rename(columns={'TF Symbol': 'TF Symbol_validated', 'TG Symbol': 'TG Symbol_validated', 'TF_type': 'TF_type_validated'})

    # Merge the non merged sentences with the found ones in candidate sentences
    merged_df = pd.merge(non_merged_df, found_sents, on=['#SentenceID', 'Sentence'], how='inner', suffixes=('_validated', '_prerenorm'))
    print(f"{len(merged_df)}/{len(non_merged_df)} rows are merged.")

    # Get non_merged columns
    non_merged_df_2 = non_merged_df[~non_merged_df['#SentenceID'].isin(merged_df['#SentenceID'])]

    # Only keep the False
    merged_df_false = merged_df[merged_df['Label_prerenorm'] == 'FALSE']

    # For False ones, we will not check NER and Norm
    merged_df_false = merged_df_false.drop(columns=[f'{T} Symbol_validated' for T in ('TF', 'TG')])

    return merged_df_false, non_merged_df_2

def fix_label_MoR(df) :
    '''
    Fix mismatch between Label and MoR in the two datasets.
    
    Diagram of all possibilities (with val=validated, pre=prerenorm):

    val_LM	pre_LM	val_V?	val_T_LM ->	V?	    T_LM
    T|A		F		T		.			F		T|A
                    F		F			T		.
                    F		T|R			F		T|R
            T|R		T		.			F		T|A
                    F		F			F		F
                            T|R			T		.
                    F		T|U			F		T|U
    F		T|A		T		.			F		F
                    F		T|A			T		.
                            T|R			F		T|R

    Summary:
        V? =  "T" if pre_LM == val_T_LM else "F"
        TM =  nan if pre_V? == "T" else (val_LM if val_T_LM is nan else val_T_LM)
    '''

    # Identify mismatches
    m_mismatch = (df['Label_prerenorm'] != df['Label_validated']) | (df['MoR_prerenorm'] != df['MoR_validated'])

    # Create combined columns for easy comparison
    df['val_LM']   = df['Label_validated'] + '|' + df['MoR_validated'].fillna('')
    df['pre_LM']   = df['Label_prerenorm'] + '|' + df['MoR_prerenorm'].fillna('')
    df['val_T_LM'] = df['true_label'].fillna('') + '|' + df['true_MoR'].fillna('')


    # Update 'Valid?' column based on mismatch condition
    df.loc[m_mismatch, 'Valid?'] = np.where(df[m_mismatch]['pre_LM'] == df[m_mismatch]['val_T_LM'], 'T', 'F')

    # Update 'T_LM' column based on validity
    df.loc[m_mismatch, 'T_LM'] = np.where(
        df.loc[m_mismatch, 'Valid?'] == 'T',  # If Valid? is 'T'
        '|',                                  # Set 'T_LM' to '|'
        np.where(                             # Otherwise 
            df.loc[m_mismatch, 'val_T_LM'] == '|',  # If 'val_T_LM' is '|'
            df.loc[m_mismatch, 'val_LM'],           # Set 'T_LM' to 'val_LM'
            df.loc[m_mismatch, 'val_T_LM']          # Otherwise, set to 'val_T_LM'
        )
    )

    # Split 'T_LM' into 'true_label' and 'true_MoR'
    df.loc[m_mismatch, ['true_label', 'true_MoR']] = df.loc[m_mismatch, 'T_LM'].str.split('|', expand=True)

    # Drop extra columns
    df.drop(columns=['Label_validated', 'MoR_validated', 'val_LM', 'pre_LM', 'val_T_LM', 'T_LM'], inplace=True)
    df.rename(columns={'Label_prerenorm': 'Label', 'MoR_prerenorm': 'MoR'}, inplace=True)

    return

def fix_NFKB_AP1_mismatches(merged_df: pd.DataFrame) -> None:
    '''
    All TF symbol mismatches should be due to NFKB/AP1. Raise assertion error if not.
    If so, remove the 'renormalisation' error (validations were done before the NFKB/AP1 fixing)
    '''

    # Get sentences validated and prerenorm sentences don't match. Use a set to ignore order in cases like 'BRAC1;BRAC2'
    tf_val_symbol_set = merged_df['TF Symbol_validated'].fillna('').str.upper().str.split(";").apply(lambda x: set(x))
    tf_pre_symbol_set = merged_df['TF Symbol_prerenorm'].fillna('').str.upper().str.split(";").apply(lambda x: set(x))
    m_mismatch = tf_val_symbol_set != tf_pre_symbol_set

    # Assert that all cases are due to NFKB / AP1 renormalisations
    assert all(merged_df[m_mismatch]['TF Symbol_prerenorm'].str.upper().isin(('NFKB', 'AP1')))
    
    # For those cases, remove 'renormalisation' error
    merged_df.loc[m_mismatch, 'TF_is_incorrect']    = np.nan
    merged_df.loc[m_mismatch, 'TF_correct_mention'] = np.nan

    # All TF symbols have been checked, so we can drop the _validated column
    merged_df.drop(columns=['TF Symbol_validated', 'TF Id_validated'], inplace=True)
    merged_df.rename(columns={'TF Symbol_prerenorm': 'TF Symbol', 'TF Id_prerenorm': 'TF Id'}, inplace=True)

    return 

def fix_TG_Symbol_ID_mismatches(merged_df, state='prerenorm'):
    '''
    Ensure TG Id/Symbol mismatches are expected
    Fix corrected normalisations
    Drop TG ID & TG Symbol columns
    '''

    # Get sentences where old and new TF Symbol don't match. Use a set to ignore order in cases like 'BRAC1;BRAC2'
    tg_pre_symbol_set = merged_df[f'TG Symbol_{state}'].fillna('').str.upper().str.split(";").apply(lambda x: set(x))
    tg_val_symbol_set = merged_df['TG Symbol_validated'].fillna('').str.upper().str.split(";").apply(lambda x: set(x))
    m_mismatch = tg_pre_symbol_set != tg_val_symbol_set

    tg_pre_id_set = merged_df[f'TG Id_{state}'].fillna('').str.split(";").apply(lambda x: set(x))
    tg_val_id_set = merged_df['TG Id_validated'].fillna('').str.split(";").apply(lambda x: set(x))
    m_id_mismatch = tg_pre_id_set != tg_val_id_set

    # Ensure the mismatch rows are the same for the symbol and the ID
    assert m_mismatch.equals(m_id_mismatch)

    # Ensure that all mismatches are due to the stage 7 of postrenormalisations
    if state == 'prerenorm':
        assert m_mismatch.sum() == 0

        # All TG symbols have been checked, so we can drop the _validated column
        merged_df.drop(columns=['TG Symbol_validated', 'TG Id_validated'], inplace=True)
        merged_df.rename(columns={'TG Symbol_prerenorm': 'TG Symbol', 'TG Id_prerenorm': 'TG Id'}, inplace=True)

    elif state == 'postrenorm':
        # Ensure mismatching TG symbols are part of the renormalised symbols
        upper_renormalised_symbols = {'CDKN1A', 'TP53', 'TRP53'}
        assert(all(merged_df[m_mismatch]['TG Symbol_postrenorm'].str.upper().isin(upper_renormalised_symbols))), "Some TG symbols have unexpected mismatches"   

        # Remove the 'normalisation' issue if present from those rows
        merged_df.loc[m_mismatch, 'TG_is_incorrect'] = np.nan
        merged_df.loc[m_mismatch, 'TG_correct_mention'] = np.nan

        # All TG symbols have been checked, so we can drop the _validated column
        merged_df.drop(columns=['TG Symbol_validated', 'TG Id_validated'], inplace=True)
        merged_df.rename(columns={'TG Symbol_postrenorm': 'TG Symbol', 'TG Id_postrenorm': 'TG Id'}, inplace=True)

    else:
        raise ValueError('Invalid state: must be postrenorm or prerenorm')
    
    return

def get_postrenorm_prerenorm_df(merged_df_valid, merged_df_false, valid_df):

    # For valid ones, we will correct the NFKB/AP1 normalisations
    fix_NFKB_AP1_mismatches(merged_df_valid)
    
    # Correct Label & MoR for both
    fix_label_MoR(merged_df_valid)
    fix_label_MoR(merged_df_false)

    # Check Ids are the same in both cases, and if so, join columns into one
    for T in ('TF', 'TG'):
        m_matching = (merged_df_false[f'{T} Id_validated'] == merged_df_false[f'{T} Id_prerenorm']) | (merged_df_false[f'{T} Id_validated'].str.contains('Complex:'))
        assert all(m_matching), f"Some {T} Ids are not the same in the validated and final df"
        merged_df_false = merged_df_false.drop(columns=f'{T} Id_prerenorm')
        merged_df_false = merged_df_false.rename(columns={f'{T} Id_validated': f'{T} Id'})

    # Get prerenorm sentences & ensure there are no TG Symbol/ID mismatches
    merged_prerenorm = merged_df_valid[~(merged_df_valid['stage'].str.contains('postnorm'))].copy()
    merged_prerenorm  = merged_prerenorm.rename(columns={'TF_type_prerenorm': 'TF_type'})    
    fix_TG_Symbol_ID_mismatches(merged_prerenorm, state='prerenorm')

    # Get postrenorm sentences & fix the TG Symbol/ID mismatches
    merged_postrenorm, _ = join_validated_df_with_valid_df(merged_df_valid.rename(columns={f'{col}_validated': col for col in ['TF_type', 'TG Id', 'TG Symbol']}), valid_df, 'postrenorm')
    fix_TG_Symbol_ID_mismatches(merged_postrenorm, state='postrenorm')
    
    # Check & drop columns
    cols = ['Label', 'MoR', 'TF Id', 'TF Symbol']
    assert all([merged_postrenorm[f'{col}_validated'].equals(merged_postrenorm[f'{col}_postrenorm']) for col in cols]), f"Some columns have mismatches"
    merged_postrenorm = merged_postrenorm.drop(columns=[f'{col}_validated' for col in cols])
    merged_postrenorm = merged_postrenorm.rename(columns={f'{col}_postrenorm': col for col in cols+['TF_type']})

    return merged_prerenorm, merged_postrenorm, merged_df_false



# PREPARE DATASETS TO VALIDATE
def select_rows_for_validation(df: pd.DataFrame, validated: pd.DataFrame = None, method: str = 'balanced', TRI_size = None, nonTRI_size = None, TFs_drawn_per_batch = 1) -> pd.DataFrame:
    '''
    Get a sample of rows to validate.

    TRI_size:    sentences with TRIs to validate
    nonTRI_size: sentences without a TRI to validate

    Method: 
        random: a randomly selected sample of TRI_size + nonTRI_size sentences
        balanced: a sample of TRI_size sentences that ensures all TFs have the same probability of being selected.

    TFs_drawn_per_batch: number of TFs to be drawn per batch in the balanced method (1 or 2)
    
    '''
    
    if method not in {'random', 'balanced'}:
        raise ValueError("Invalid method. Choose 'random' or 'balanced'.")
    

    if method == 'random':
        # Method 1: Random Sampling (non-TRI & TRI sentences)
        if validated is not None:
            m =~df['#SentenceID'].isin(validated['#SentenceID'])
            df = df[m]
        
        TRI_sample = df[df['Valid'] == 'Valid'].sample(TRI_size)
        nonTRI_sample = df[df['Valid'] == 'Non valid'].sample(nonTRI_size)
        to_validate = pd.concat([TRI_sample, nonTRI_sample])
        print(f'We will validate {len(to_validate)} sentences with the {method} method')


    elif method == 'balanced':
        # Method 2: Balanced Sampling based on TF Id frequency (only TRI sentences)
        TRI_df = df[df['Valid'] == 'Valid']
        sorted_TFs = TRI_df['TF Id'].value_counts().index
        batch_size = math.ceil(len(sorted_TFs) / TRI_size) * TFs_drawn_per_batch
        batched_TFs = [sorted_TFs[i:i + batch_size] for i in range(0, len(sorted_TFs), batch_size)]
        to_validate = pd.DataFrame()

        if validated is not None:
            m = (validated['Valid'] == 'Valid') if 'Valid' in validated.columns else (validated['Label'] == 'TRUE')
            validated_TFs_set = set(validated[m]['TF Id'])
        else:
            validated_TFs_set = set()

        for batch in batched_TFs:
            validated_count = sum(TF in validated_TFs_set for TF in batch)

            if validated_count < TFs_drawn_per_batch:
                # Get 1 random sentence per TF ID in the batch
                m = TRI_df['TF Id'].isin(batch)
                filtered_df = TRI_df[m].sample(frac=1).drop_duplicates(subset='TF Id')

                # If one TF is already validated, remove it and only draw 1 sentence
                if validated_count == 1:
                    # Remove already validated TFs from the filtered_df
                    filtered_df = filtered_df[~filtered_df['TF Id'].isin(validated_TFs_set)]
                    sample_size = min(1, len(filtered_df))  # Draw only 1 sentence in this case
                else:
                    # Draw the minimum of TFs_drawn_per_batch or available filtered_df rows
                    sample_size = min(TFs_drawn_per_batch, len(filtered_df))

                # Add the appropriate number of sentences to validate
                if sample_size > 0:
                    to_validate = pd.concat([to_validate, filtered_df.sample(sample_size, replace=False)])

        print(f'We will validate {len(to_validate)}/{len(batched_TFs*TFs_drawn_per_batch)} sentences with the {method} method, with {TFs_drawn_per_batch} TFs selected per batch')

    # Specify the method used in the validated sentences
    to_validate['method'] = method

    # Prepare for Excel
    reorder_df_for_Excel(to_validate, TaxID_to_species)

    return to_validate

def reorder_df_for_Excel(to_validate: pd.DataFrame, TaxID_to_species: dict = TaxID_to_species) -> pd.DataFrame:
    '''Reorder the columns of the dataframe in the way to be shown in the Excel used for manual validation'''

    # Add span_sentence by <TF>...</TF>
    to_validate['span_sentence'] = to_validate.apply(lambda row: row['Sentence'].replace('[TF]', f"<TF>{row['TF']}</TF>").replace('[TG]', f"<TG>{row['TG']}</TG>"), axis=1)

    # Rename columns & cell names
    to_validate = to_validate.rename(columns={"Valid": "Label"})
    to_validate.loc[to_validate['Label'] == 'Valid', 'Label'] = 'TRUE'
    to_validate.loc[to_validate['Label'] == 'Non valid', 'Label'] = 'FALSE'
    to_validate.loc[to_validate['Label'] == 'FALSE', 'MoR'] = ''

    # Add extra columns
    extra_cols = ['Valid?', 'Incorrect? TF', 'Incorrect? TG', 'Correct name TF', 'Correct name TG', 'true_label', 'true_MoR', 'dir-*', 'Other issues', 'OI2', 'Explanation', ]
    for extra_col in extra_cols:
        to_validate[extra_col] = ''

    # Translate TaxID into species
    def translate_taxid_to_species(TaxID):
        return ';'.join(TaxID_to_species.get(id, 'other') for id in TaxID.split(';'))
    to_validate['TF Species'] = to_validate['TF TaxID'].apply(translate_taxid_to_species)
    to_validate['TG Species'] = to_validate['TG TaxID'].apply(translate_taxid_to_species)

    ordered_main_cols = ['#SentenceID', 'TF', 'TG', 'span_sentence', 
        'TF Symbol', 'Incorrect? TF', 'Correct name TF', 
        'TG Symbol', 'Incorrect? TG', 'Correct name TG', 
        'Label', 'MoR', 'Valid?', 'dir-*', 'true_label', 'true_MoR', 'Other issues', 'OI2', 'Explanation', 
        'TF Species', 'TG Species',
        'Valid score', 'MoR scores']
    remaining_cols = to_validate.columns.difference(ordered_main_cols)
    to_validate = to_validate[ordered_main_cols + list(remaining_cols)]

    return to_validate


# Show validated datasets
def get_PMIDs_used_in_training(training_path: str) -> set:
    '''Return a set of the PMIDs used for training the classifiers (to exclude from validation)'''
    training_data = pd.read_csv(training_path, sep='\t', header=0)
    training_data['PMID'] = training_data['#TRI ID'].apply(lambda row: row.split(':')[0])
    training_PMIDs = set(training_data['PMID'])
    print(f"Training PMIDs:\t We used {len(training_PMIDs)} PMIDs in the training set, that we will filter out for validation.")
    return training_PMIDs

def load_valid_nonvalid_df(valid_path: str, nonvalid_path: str, training_path: str) -> tuple:
    '''Return valid and nonvalid(subset) dfs, with the PMIDs used in training filtered out'''
    # Load and concatenate
    valid_df    = pd.read_csv(valid_path, sep='\t', header=0, dtype=str)
    nonvalid_df = pd.read_csv(nonvalid_path, sep='\t', header=0, dtype=str)

    # Filter PMIDs used in training
    training_PMIDs = get_PMIDs_used_in_training(training_path)
    valid_df    = valid_df[~valid_df['PMID'].isin(training_PMIDs)]
    nonvalid_df = nonvalid_df[~nonvalid_df['PMID'].isin(training_PMIDs)]

    # Remove hash from #SentenceID
    valid_df['#SentenceID']     = valid_df['#SentenceID'].apply(lambda row: ":".join(row.split(':')[i] for i in [1,4,5,6]))
    nonvalid_df['#SentenceID']  = nonvalid_df['#SentenceID'].apply(lambda row: ":".join(row.split(':')[i] for i in [1,4,5,6]))

    return valid_df, nonvalid_df 

def display_validated_per_stage(validated_df: pd.DataFrame) -> str: 
    '''Show in markdown format the number of validated sentences in each stage'''

    # Get length
    s1 = (validated_df['stage'] == 's1').sum()
    s2 = (validated_df['stage'] == 's2').sum()
    s3 = (validated_df['stage'] == 's3').sum()
    s4_db_b = (validated_df['stage'] == 's4_dbTF').sum()


    s4_co_b = ((validated_df['stage'] == 's4_coTF') & (validated_df['method'] == 'balanced')).sum()
    s4_co_r = ((validated_df['stage'] == 's4_coTF') & (validated_df['method'] == 'random')).sum()
    s5_co_b = ((validated_df['stage'] == 's5_coTF') & (validated_df['method'] == 'balanced')).sum()
    s5_co_r = ((validated_df['stage'] == 's5_coTF') & (validated_df['method'] == 'random')).sum()
    s6_co_b = ((validated_df['stage'] == 's6_coTF') & (validated_df['method'] == 'balanced')).sum()
    s6_co_r = ((validated_df['stage'] == 's6_coTF') & (validated_df['method'] == 'random')).sum()

    s6_co_ll_b = ((validated_df['stage'] == 's6_coTF_ll') & (validated_df['method'] == 'balanced')).sum()
    s6_co_ll_r = ((validated_df['stage'] == 's6_coTF_ll') & (validated_df['method'] == 'random')).sum()

    # Get type
    s1t = ",".join(validated_df[validated_df['stage'] == 's1']['TF_type'].unique())
    s2t = ",".join(validated_df[validated_df['stage'] == 's2']['TF_type'].unique())
    s3t = ",".join(validated_df[validated_df['stage'] == 's3']['TF_type'].unique())

    s1m = (validated_df[validated_df['stage'] == 's1']['Mutated_TF'] == 'TRUE').sum()
    s1f = (validated_df[validated_df['stage'] == 's1']['Label'] == 'FALSE').sum()
    s2f = (validated_df[validated_df['stage'] == 's2']['Label'] == 'FALSE').sum()
    s3f = (validated_df[validated_df['stage'] == 's3']['Label'] == 'FALSE').sum()

    vdbTF = (validated_df['TF_type'] == 'dbTF').sum()
    vGTF = (validated_df['TF_type'] == 'GTF').sum()
    vcoTF = (validated_df['TF_type'] == 'coTF').sum()
    vcoTF_ll = (validated_df['TF_type'] == 'coTF_ll').sum()

    md_text = (f'''
    * Proof of concept run on 10% of PubMed
        * **Stage 1**: {s1} randomly selected {s1t} sentences, {s1m} of which mutated, and {s1f} false.
    * First full run on all PubMed
        * **Stage 2:** {s2} randomly selected {s2t} sentences, {s2f} false.
        * **Stage 3:** {s3} randomly selected {s3t} sentences, {s3f} false.
    * Second full run on all PubMed
        * **Stage 4**:
            * {s4_db_b} balanced dbTF sentences, 
            * {s4_co_r} random coTF sentences,
            * {s4_co_b} balanced coTF sentences. 
        * **Stage 5**:
            * {s5_co_r} random coTF sentences,
            * {s5_co_b} balanced coTF sentences. 
        * **Stage 6**:
            * {s6_co_r} random coTF sentences
            * {s6_co_b} balanced coTF sentences 
            * {s6_co_ll_r} random less likely coTF sentences
            * {s6_co_ll_b} balanced less likely coTF sentences

    In total, we have {vdbTF} dbTF, {vGTF} GTF, {vcoTF} coTF  and {vcoTF_ll} less likely coTF validated sentences.
    ''')

    return md_text.replace('\n', '<br>')
