import numpy as np
from score.density_estimate import DensityEstimate
from features.features import Features
from glob import glob
import os
import pandas as pd

def read_stats(stats_root, features):
    stats = {}
    for dist in ['native', 'reference']:
        for interaction in features:
            fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
            assert os.path.exists(fname), fname
            if interaction not in stats: stats[interaction] = {}
            stats[interaction][dist] = DensityEstimate.read(fname)
    return stats

def pair_features(protein, data_root, pairs_root):
    interactions = ['hbond',  'saltbridge', 'contact', 'shape', 'mcss']
    features = Features(data_root + '/' + protein, max_poses=100)
    features.load_features(interactions)
    features = features.get_view()

    # get cross-docked ligands
    ligands = []
    for ligand in sorted(features['gscore'].keys()):
        if 'native' in ligand: continue
        lig, grid = ligand.split('-to-')
        if '_lig' in lig:
            lig = lig.replace('_lig', '')
        if lig != grid and 'CHEMBL' not in lig:
            ligands += [ligand]

    df = []
    for i, ligand1 in enumerate(ligands):
        for ligand2 in ligands[i+1:]:
            for r1 in range(len(features['gscore'][ligand1])):
                for r2 in range(len(features['gscore'][ligand2])):
                    feats = [features[interaction][(ligand1, ligand2)][r1, r2]
                             for interaction in interactions]
                    gscore1 = features['gscore'][ligand1][r1]
                    gscore2 = features['gscore'][ligand2][r2]
                    rmsd1 = features['rmsd'][ligand1][r1]
                    rmsd2 = features['rmsd'][ligand2][r2]
                    df += [[protein,
                            ligand1, ligand2,
                            r1, r2,
                            gscore1, gscore2,
                            rmsd1, rmsd2]
                            + feats]
    df =  pd.DataFrame(df, columns=['protein',
                                    'ligand1', 'ligand2',
                                    'rank1', 'rank2',
                                    'gscore1', 'gscore2',
                                    'rmsd1', 'rmsd2']
                                     +interactions)
    df.to_csv('{}/{}.csv'.format(pairs_root, protein), index=False)

def compute_stats(protein, pairs_root, stats_root, features):
    df = pd.read_csv('{}/{}.csv'.format(pairs_root, protein))
    for feature in features:
        if feature == 'mcss':
            sd = 0.03*6
            domain = (0, 6)
        else:
            sd = 0.03
            domain = (0, 1)
        
        nat_vals = df.loc[(df.rmsd1 <= 2.0)&(df.rmsd2 <= 2.0), feature]
        ref_vals = df.loc[:, feature]
        nat = DensityEstimate(domain=domain, sd=sd).fit(nat_vals)
        ref = DensityEstimate(domain=domain, sd=sd).fit(ref_vals)
        nat.write('{}/{}/native_{}.de'.format(stats_root, protein, feature))
        ref.write('{}/{}/reference_{}.de'.format(stats_root, protein, feature))

def merge_stats(proteins, stats_root, merged_stats_fname, features):
    for feature in features:
        nat_des, ref_des = [], []
        for protein in proteins:
            nat_fname = '{}/{}/native_{}.de'.format(stats_root, protein, feature)
            ref_fname = '{}/{}/reference_{}.de'.format(stats_root, protein, feature)
            nat_des += [DensityEstimate.read(nat_fname)]
            ref_des += [DensityEstimate.read(ref_fname)]

        nat_fname = merged_stats_fname.format('native', feature)
        ref_fname = merged_stats_fname.format('reference', feature)
        DensityEstimate.merge(nat_des).write(nat_fname)
        DensityEstimate.merge(ref_des).write(ref_fname)
