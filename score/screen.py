import numpy as np
import pandas as pd
from utils import np_load
from schrodinger.structure import StructureReader, StructureWriter

def load_features_screen(features, root):
    single = np.load(f'{root}/gscore1.npy')

    raw = {}
    for feature in features:
        raw[feature] = np_load(f'{root}/{feature}.npy')
    return single, raw

def scores_to_csv(pv, out):
    """
    Write docking and ComBind scores to text.
    """
    titles, glide, combind = [], [], []
    with StructureReader(pv) as reader:
        next(reader)
        for st in reader:
            titles += [st.title]
            glide += [st.property['r_i_docking_score']]
            combind += [st.property['r_i_combind_score']]

    df = pd.DataFrame(np.vstack([titles, glide, combind]).T,
                      columns = ['ID', 'GLIDE', 'COMBIND'])
    df.to_csv(out, index=False)

def apply_scores(pv, scores, out):
    """
    Add ComBind screening scores to a poseviewer.
    """

    scores = np.load(scores)

    with StructureReader(pv) as reader, StructureWriter(out) as writer:
        st = next(reader)
        st.property['r_i_combind_score'] = 1000.0
        writer.append(st)
        for st, score in zip(reader, scores):
            st.property['r_i_combind_score'] = score
            writer.append(st)

def screen(single, raw, stats, alpha):
    energies = {}
    for feature in raw:
        _raw = raw[feature]
        _stats = stats[feature]
        energies[feature] = (  np.log(_stats['native'](_raw))
                             - np.log(_stats['reference'](_raw)))

    pair_energy = 0
    for feature, energy in energies.items():
        pair_energy += energy

    n = pair_energy.shape[1]

    alpha /= 0.5 * n / (1 + (n-1)*0.5)

    pair_energy = pair_energy.mean(axis=1)
    combind_energy = pair_energy/alpha - single
    return combind_energy
