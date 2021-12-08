import sqlite3
import os
import pandas as pd
import numpy as np
from schrodinger.structure import SmilesStructure
import click

# Macrocycles are hard to dock.
MACROCYCLE_THRESH = 8

class CHEMBLDB:
    def __init__(self, chembldb, uniprot_chembl):
        # CHEMBL database was downloaded from the chembl website in sql format.
        self.conn = sqlite3.connect('file:{}?mode=ro'.format(chembldb), uri=True)
        self.cur = self.conn.cursor()
        
        # File of links from uniprot id to chembl id was downloaded from uniprot.
        # I just searched for all uniprot entriies with a chembl cross-reference
        # available.
        self.uniprot_chembl = pd.read_csv(uniprot_chembl, sep='\t', index_col=0)
        self.uniprot_chembl = self.uniprot_chembl['Cross-reference (ChEMBL)'].apply(lambda x: x.strip(';').split(';'))

    def __enter__ (self):
        return self
    def __exit__ (self, *exc):
        self.conn.close()

    def chembl_to_tid(self, chembl):
        self.cur.execute("SELECT tid FROM target_dictionary WHERE chembl_id=?", (chembl,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def tid_to_target_type(self, tid):
        self.cur.execute("SELECT target_type FROM target_dictionary WHERE tid=?", (tid,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def tid_to_assays(self, tid, protein_complex, homologous):
        if protein_complex and homologous:
            confidence = '(confidence_score=6 OR confidence_score=7)'
        elif protein_complex and not homologous:
        	confidence = 'confidence_score=7'
        elif not protein_complex and homologous:
        	confidence = '(confidence_score=8 OR confidence_score=9)'
        else:
        	confidence = 'confidence_score=9'

        self.cur.execute("SELECT assay_id FROM assays WHERE tid=? AND "+confidence, (tid,))
        return [row[0] for row in self.cur.fetchall()]

    def assay_to_chemblid(self, assay):
        self.cur.execute("SELECT chembl_id FROM assays WHERE assay_id=?", (assay,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def assay_to_molregnos(self, assay):
        self.cur.execute("SELECT molregno FROM activities WHERE assay_id=?", (assay,))
        return [row[0] for row in self.cur.fetchall()]

    def molregno_to_smiles(self, molregno):
        self.cur.execute("SELECT canonical_smiles FROM compound_structures WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def molregno_to_molw(self, molregno):
        self.cur.execute("SELECT mw_freebase FROM compound_properties WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        if not len(rows):
            return 0
        return rows[0][0]

    def molregno_to_smiles(self, molregno):
        self.cur.execute("SELECT canonical_smiles FROM compound_structures WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        if not rows:
            return None
        return rows[0][0]

    def molregno_to_chemblid(self, molregno):
        self.cur.execute("SELECT chembl_id FROM molecule_dictionary WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def molregno_and_assay_to_activities(self, molregno, assay):
        self.cur.execute("SELECT standard_type, standard_value, standard_units, relation, activity_comment FROM activities WHERE molregno=? AND assay_id=?", (molregno, assay))
        return self.cur.fetchall()

    def chembl_to_activities(self, chembl, protein_complex, homologous):
        # chemblID, SMILES, MOLW, affinity
        activities = []
        tid = self.chembl_to_tid(chembl)
        for assay in self.tid_to_assays(tid, protein_complex, homologous):
            assay_chembl_id = self.assay_to_chemblid(assay)
            for molregno in self.assay_to_molregnos(assay):
                molw = self.molregno_to_molw(molregno)
                smiles = self.molregno_to_smiles(molregno)
                chembl_id = self.molregno_to_chemblid(molregno)
                for activity in self.molregno_and_assay_to_activities(molregno, assay):
                    activities += [[assay_chembl_id, chembl_id, molw, smiles] + list(activity)]
        return pd.DataFrame(activities,
                            columns=['assay_chembl_id', 'ligand_chembl_id', 'mw_freebase',
                                     'canonical_smiles', 'standard_type', 'standard_value',
                                     'standard_units', 'relation', 'comment'])

    def uniprot_to_chembl(self, uniprot):
        if uniprot not in self.uniprot_chembl:
            return None

        for chembl_id in self.uniprot_chembl.loc[uniprot]:
            tid = self.chembl_to_tid(chembl_id)
            target_type = self.tid_to_target_type(tid)
            if target_type == 'SINGLE PROTEIN':
                return chembl_id

################################################################################

def get_chembl_id(uniprot, chembldb, uniprot_chembl):
    with CHEMBLDB(chembldb, uniprot_chembl) as chembldb:
        chembl = chembldb.uniprot_to_chembl(uniprot)
    return chembl

def standardize_nonbinders(activities, affinity_thresh):
    # Set 'Not Active's to affinity_thresh
    activities.loc[activities['comment'].isin([None]), 'comment'] = ''
    duds = [('Not Active' in s) for s in activities['comment']]
    duds = np.array(duds)
    activities.loc[duds, 'standard_units'] = 'nM'
    activities.loc[duds, 'standard_value'] = affinity_thresh
    activities.loc[duds, 'relation'] = '='

    # Most nonbinders don't have equality relation.
    mask  = activities['standard_value'] >= affinity_thresh
    mask *= activities['relation'].isin(['>', '>='])
    activities.loc[mask, 'relation'] = '='

    # Cap affinity values.
    mask  = activities['standard_value'] >= affinity_thresh
    activities.loc[mask, 'standard_value'] = affinity_thresh
    return activities

def get_activities(chembl, chembldb, uniprot_chembl, protein_complex, homologous, affinity_thresh):
    with CHEMBLDB(chembldb, uniprot_chembl) as chembldb:
        activities = chembldb.chembl_to_activities(chembl, protein_complex, homologous)
    activities['target_chembl_id'] = chembl

    # Standardize units to nM
    m = {'M': 10**9, 'mM': 10**6, 'uM': 10**3, 'pM': 10**-3}
    for unit, relation in m.items():
        mask = activities['standard_units'] == unit
        activities.loc[mask, 'standard_value'] *= relation
        activities.loc[mask, 'standard_units'] = 'nM'

    activities = standardize_nonbinders(activities, affinity_thresh)

    return activities

################################################################################

def filter_activities(activities, activity_type, molw_thresh):
    if activity_type == 'all':
        activity_types = ['IC50', 'Ki', 'Kd']
    else:
        activity_types = [activity_type]
    mask = activities['standard_type'].isin(activity_types)
    print('Removing {} rows b/c standard_type not in {}'.format(len(mask)-sum(mask),
                                                                 activity_types))
    print('Set of offending values is {}'.format(set(activities[~mask]['standard_type'])))
    activities = activities.loc[mask]

    mask = activities['standard_value'].notna()
    print('Removing {} rows b/c standard_value is na'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['standard_value'] != 0
    print('Removing {} rows b/c standard_value is 0'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = ~activities['canonical_smiles'].isin([None])
    print('Removing {} rows b/c canonical_smiles is None'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['mw_freebase'] <= molw_thresh+100 # Not desalted yet.
    print('Removing {} rows b/c mw_freebase > {}'.format(len(mask)-sum(mask),
                                                         molw_thresh+100))
    activities = activities.loc[mask]

    mask = activities['standard_units'] == 'nM'
    print('Removing {} rows b/c standard_units != nM'.format(len(mask)-sum(mask)))
    print('Set of offending values is {}'.format(set(activities[~mask]['standard_units'])))
    activities = activities.loc[mask]

    mask = activities['relation'] == '='
    print('Removing {} rows b/c relation != ='.format(len(mask)-sum(mask)))
    print('Set of offending values is {}'.format(set(activities[~mask]['relation'])))
    activities = activities.loc[mask]

    return activities

################################################################################

def desalt(smiles):
    ligand = ''
    for molecule in smiles.split('.'):
        if len(molecule) > len(ligand):
            ligand = molecule
    return ligand

def get_structure(smiles):
    smi = SmilesStructure(smiles)
    try:
        st = smi.get3dStructure(True)
        stereo = True
    except:
        try:
            st = smi.get3dStructure(False)
            stereo = False
        except:
            print('Error processing {}'.format(smiles))
            return SmilesStructure('C').get3dStructure(), False

    return st, stereo

def is_macrocycle(st):
    ring_sizes = [0]+[len(ring.atom) for ring in st.ring]
    return max(ring_sizes) > MACROCYCLE_THRESH

def _get_properties(smiles):
    properties = {}
    properties['SMILES'] = desalt(smiles)
    st, properties['STEREO'] = get_structure(properties['SMILES'])
    properties['MACROCYCLE'] = is_macrocycle(st)
    properties['MOLW'] = st.total_weight
    return pd.Series(properties)

def get_properties(activities):
    properties = activities.canonical_smiles.apply(_get_properties)
    return pd.concat([activities, properties], axis=1)

def filter_properties(activities, ambiguous_stereo, molw_thresh):
    mask = activities['MOLW'] <= molw_thresh
    print('Removing {} rows b/c molw > {}'.format(len(mask)-sum(mask),
                                                  molw_thresh))
    activities = activities.loc[mask]

    mask = ~activities['MACROCYCLE']
    print('Removing {} rows b/c macrocycle'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    if not ambiguous_stereo:
    	mask = activities['STEREO']
    	print('Removing {} rows b/c ambiguous stereochemistry'.format(len(mask)-sum(mask)))
    	activities = activities.loc[mask]
    return activities

################################################################################

def collapse_duplicates(activities, seperate_activity_types):
    if seperate_activity_types:
        keys = ['SMILES', 'standard_type']
    else:
        keys = ['SMILES']

    averages = activities.loc[:, keys+['standard_value']].groupby(keys).mean()
    activities = activities.groupby(keys, as_index=False).first()
    activities['standard_value'] = [averages.loc[tuple([row[key] for key in keys])]['standard_value']
                                    for _, row in activities.iterrows()]
    return activities

################################################################################

@click.group()
def main():
    pass

@main.command()
@click.option('--protein-complex', is_flag=True)
@click.option('--homologous', is_flag=True)
@click.option('--ambiguous-stereo', is_flag=True)
@click.option('--activity-type', default='all')
@click.option('--affinity-thresh', default=10000)
@click.option('--molw-thresh', default=500)
@click.option('--output-fname')
@click.argument('uniprot_or_chembl')
@click.argument('chembldb', default='/oak/stanford/groups/rondror/users/jpaggi/pldb_data/raw/chembl_27.db')
@click.argument('uniprot_chembl', default='/oak/stanford/groups/rondror/users/jpaggi/pldb_data/raw/uniprot-chembl.tsv')
def query(protein_complex, homologous, ambiguous_stereo, activity_type,
          affinity_thresh, molw_thresh, output_fname,
          uniprot_or_chembl, chembldb, uniprot_chembl):
    """
    protein_complex (bool): protein complexes will necessarily have a lower
        confidence score in CHEMBL, so when looking for data for a complex
        you have to lower the threshold.
    homologous (bool): allow matches with a confidence score of 8. This means
        that they weren't able to map the activity to a particular version of a
        protein. For instance, if a paper reports "binding affinity to the D2
        dopamine receptor", but doesn't specify human, rat, etc..
    ambiguous_stereo (bool): allow ligands where the stereochemistry isn't
        completely deteremined.
    activity_type (str: IC50, Ki, Kd, EC50): specifies what kinds of measurements
        to accept. If not specified, will give [IC50, Ki, Kd].
    affinity_thresh (float): affinity value to assign to "Not actives" and
        minimum affinity for which to accept ligands with a "greater than"
        relation. 10000 (10 micro-molar) is generally a good option because this
        is a standard threshold for binding in screens.
    molw-thresh (float): specifies how large of molecules to include in the
        results. 500 Da is a reasonable threshold to get "small molecules"
    output_fname (str): where to put the results. By default,
        CHEMBLID_ACTIVITYTYPE.csv
    uniprot_or_chembl (str): target for which to get affinities.
    chembldb (str): path to the chembl database
    uniprot_chembl (str): path to file cross-referencing uniprot to chembl.
    """

    # Resolve chembl target id.
    if uniprot_or_chembl[:6] == 'CHEMBL':
        chembl = uniprot_or_chembl
    else:
        chembl = get_chembl_id(uniprot_or_chembl, chembldb, uniprot_chembl)

    if chembl is None:
        print('No CHEMBL entry...')
        if output_fname:
            with open(output_fname, 'w') as fp:
                fp.write('No CHEMBL entry...\n')
        return

    if output_fname is None:
        output_fname = '{}_{}.csv'.format(chembl, activity_type)

    # Look up affinities.
    activities = get_activities(chembl, chembldb, uniprot_chembl,
                                protein_complex, homologous, affinity_thresh)
    if activities.shape[0]:
        activities = filter_activities(activities, activity_type, molw_thresh)

    if activities.shape[0]:
        activities = get_properties(activities)
        activities = filter_properties(activities, ambiguous_stereo, molw_thresh)

    if not activities.shape[0]:
        print('No valid ligands...')
        if output_fname:
            with open(output_fname, 'w') as fp:
                fp.write('No valid ligands...\n')
        return

    activities = activities.sort_values(['assay_chembl_id', 'ligand_chembl_id'])
    activities.to_csv(output_fname, index=False)

@main.command()
@click.argument('input_csv')
@click.argument('output_csv')
@click.option('--seperate-activity-types', is_flag=True)
@click.option('--IC50-factor', default=2.3)
def unique(input_csv, output_csv, seperate_activity_types, ic50_factor):
    """
    If seperate_activity_types, don't merge e.g. IC50's and Ki's.
    """
    activities = pd.read_csv(input_csv)
    
    mask = (activities.standard_type == 'IC50') & (activities.standard_value >= 1000)
    activities.loc[mask, 'standard_value'] /= ic50_factor

    # Since we're going to merge all the values, assay_chembl_id will no longer
    # make sense.
    activities = activities.drop(columns='assay_chembl_id')
    activities = collapse_duplicates(activities, seperate_activity_types)
    activities.to_csv(output_csv, index=False)

if __name__ == '__main__':
    main()
