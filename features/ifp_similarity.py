import pandas as pd
import numpy as np

def merge_hbonds(ifp):
    """
    Reads IFP file and merges hbond acceptors and donors.

    Setting the label to hbond for the hbond_donors and hbond_acceptors while
    changing the residue names allows for only donor+donor or acceptor+acceptor
    to be counted as overlapping, but them to be merged into the same similarity
    measure.
    """

    mask = ifp.label=='hbond_acceptor'
    ifp.loc[mask, 'protein_res'] = [res+'acceptor' for res in ifp.loc[mask, 'protein_res']]
    ifp.loc[mask, 'label'] = 'hbond'
    
    mask = ifp.label=='hbond_donor'
    ifp.loc[mask, 'protein_res'] = [res+'donor' for res in ifp.loc[mask, 'protein_res']]
    ifp.loc[mask, 'label'] = 'hbond'
    return ifp

def ifp_tanimoto(ifps1, ifps2, feature):
    """
    Computes the tanimoto distance between ifp1 and ifp2 for feature.
    """
    if feature == 'hbond':
        ifps1 = [merge_hbonds(ifp) for ifp in ifps1]
        ifps2 = [merge_hbonds(ifp) for ifp in ifps2]

    ifps1 = [ifp.loc[ifp.label == feature] for ifp in ifps1]
    ifps2 = [ifp.loc[ifp.label == feature] for ifp in ifps2]

    ifps1 = [ifp.set_index('protein_res') for ifp in ifps1]
    ifps2 = [ifp.set_index('protein_res') for ifp in ifps2]

    sims = np.zeros((len(ifps1), len(ifps2)))
    for i, ifp1 in enumerate(ifps1):
        for j, ifp2 in enumerate(ifps2):
            total = ifp1['score'].sum() + ifp2['score'].sum()
            overlap = ifp1.join(ifp2, rsuffix='_2', how='inner')
            overlap = overlap['score']**0.5 * overlap['score_2']**0.5
            overlap = overlap.sum()

            sims[i, j] = (1 + overlap) / (2 + total - overlap)
    return sims
