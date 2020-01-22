import numpy as np
import pandas as pd

def rev_comp(seq):
    """
    Gets the reverse complement of a sequence
    """
    match_dict = {'A': 'T',
                  'T': 'A',
                  'C': 'G',
                  'G': 'C'}

    return ''.join([match_dict[x] for x in seq][::-1])

OUTFILE = ''
FULL_LENGTH_FILE = ''
CLEAVAGE_STRAND_FILE = ''

# import count data
counts = pd.read_csv(FULL_LENGTH_FILE, delim_whitespace=True, header=None)
counts.columns = ['sequence','count']
counts = counts.set_index('sequence')

counts_short = pd.read_csv(CLEAVAGE_STRAND_FILE, delim_whitespace=True, header=None)
counts_short.columns = ['sequence','count']
counts_short['len'] = [len(x) for x in counts_short['sequence']]
counts_short = counts_short.set_index('sequence')
counts5p = counts_short[counts_short['len'] < 12]
counts3p = counts_short[counts_short['len'] >= 12]

# get version of 3p counts that truncates the first 2 nucleotides
counts3p_trunc = counts3p.copy()
counts3p_trunc['truncated'] = [x[:-2] for x in counts3p.index]
print(len(counts3p_trunc))
counts3p_trunc = counts3p_trunc.drop_duplicates(subset='truncated')
print(len(counts3p_trunc))
counts3p_trunc = counts3p_trunc.reset_index().set_index('truncated')

# iterates through sequences and records matching sequences into a dataframe
goal_df = []
all_sequences = list(counts.index)
for ix, guide in enumerate(all_sequences):
    if (ix % 100) == 0:
        print(ix)
    guide_rc = rev_comp(guide)
    passengers = []
    for passenger in all_sequences[ix:]:
        if guide_rc[2:] == passenger[:-2]:
            passengers.append(passenger)
    if len(passengers) > 1:
        continue
        
    elif len(passengers) == 1:
        passenger = passengers[0]
    
    else:
        passenger = None
    
    if guide[:-12] in counts5p.index:
        guide5p = guide[:-12]
    else:
        guide5p = None 
    
    if guide[-12:] in counts3p.index:
        guide3p = guide[-12:]
    else:
        guide3p = None
    
    if passenger is not None:
        if passenger[:-12] in counts5p.index:
            pass5p = passenger[:-12]
        else:
            pass5p = None 

        if passenger[-12:] in counts3p.index:
            pass3p = passenger[-12:]
        else:
            pass3p = None
    else:
        if guide_rc[2:-10] in counts5p.index:
            pass5p = guide_rc[2:-10]
        else:
            pass5p = None
        
        if guide_rc[-10:] in counts3p_trunc.index:
            pass3p = counts3p_trunc.loc[guide_rc[-10:]]['sequence']
        else:
            pass3p = None
        
    goal_df.append([guide, passenger, guide5p, guide3p, pass5p, pass3p])

# convert to dataframe and add column names
goal_df = pd.DataFrame(goal_df)
goal_df.columns = ['guide','passenger','guide5p','guide3p','pass5p','pass3p']

# add count information back in
goal_df_with_counts = goal_df.copy()
for c, df in zip(goal_df.columns, [counts, counts, counts5p, counts3p, counts5p, counts3p]):
    goal_df_with_counts['{}_count'.format(c)] = [df.loc[x]['count'] if x is not None else None for x in goal_df[c]]

# write to outfile
goal_df_with_counts.to_csv(OUTFILE, sep='\t', index=False)
