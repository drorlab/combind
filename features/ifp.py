"""
Compute interaction fingerprints for poseviewer files.
"""

import tempfile
import os
import click
import numpy as np
import pandas as pd
from rdkit.Chem import MolFromSmarts
from rdkit.Chem.rdmolfiles import MaeMolSupplier
import gzip

################################################################################

def _convert_mae(original_mae, converted_mae, poses):
    from schrodinger.structure import StructureReader, StructureWriter
    with StructureReader(original_mae) as sts, StructureWriter(converted_mae) as writer:
        for i, st in enumerate(sts):
            if i > poses:
                break
            for k in st.property.keys():
                if 'title' not in k:
                    st.property.pop(k)
            writer.append(st)

def convert_mae(original_mae, converted_mae, poses=float('inf')):
    imp = 'import sys; sys.path.append("{}"); import ifp'.format(os.path.dirname(os.path.abspath(__file__)))
    cmd = 'ifp._convert_mae("{}", "{}", {})'.format(original_mae, converted_mae, poses)
    os.system('run python3 -c \'{};{}\''.format(imp, cmd))

################################################################################
def resname(atom):
    info = atom.GetPDBResidueInfo()
    if info is None:
        return ''
    return ':'.join(map(lambda x: str(x).strip(),
                        [info.GetChainId(), str(info.GetResidueNumber()),
                         info.GetResidueName(), info.GetInsertionCode()]))

def atomname(atom):
    pdb = atom.GetPDBResidueInfo()
    if pdb is None:
        return str(atom.GetIdx())
    return pdb.GetName().strip()

def coords(atom):
    return atom.GetOwningMol().GetConformer(0).GetAtomPosition(atom.GetIdx())

def centroid_coords(atoms):
    _coords = np.array([coords(atom) for atom in atoms])
    _coords = _coords.mean(axis=0)
    return _coords

def distance(atom1, atom2):
    return coords(atom1).Distance(coords(atom2))

def angle_atom(atom1, atom2, atom3):
    v1 = coords(atom1) - coords(atom2)
    v3 = coords(atom3) - coords(atom2)
    return v1.AngleTo(v3) * 180.0 / np.pi

def angle_vector(v1, v2):
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    angle =  np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    angle *= 180 / np.pi
    if angle > 90:
        angle = 180 - angle
    assert 0 <= angle <= 90, angle
    return angle

class Molecule:
    def __init__(self, mol, is_protein, settings):
        self.mol = mol
        self.is_protein = is_protein
        self.settings = settings

        self.pipi = self.init_pipi()
        self.contacts = self.init_contacts()
        self.hbond_donors, self.hbond_acceptors = self.init_hbond()
        self.charged, self.charge_groups = self.init_saltbridge()

    def init_contacts(self):
        coord, vdw, atom_name, res_name = [], [], [], []
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() not in self.settings['nonpolar']: continue
            coord += [coords(atom)]
            atom_name += [atomname(atom)]
            res_name += [resname(atom)]
            vdw += [self.settings['nonpolar'][atom.GetAtomicNum()]]
        
        if coord:
            coord = np.vstack(coord)
            vdw = np.array(vdw)
        return coord, vdw, res_name, atom_name

    def init_pipi(self):
        rings = self.get_aromatic_rings()
        centroids, normals, atom_name, res_name = [], [], [], []
        for ring in rings:
            centroids += [self.get_centroid(ring)]
            normals += [self.get_normal(ring)]
            res_name += [resname(self.mol.GetAtomWithIdx(ring[0]))]
            atom_name += [','.join([atomname(self.mol.GetAtomWithIdx(r)) for r in ring])]
        
        if centroids:
            centroids = np.vstack(centroids)
            normals = np.vstack(normals)
        return centroids, normals, res_name, atom_name

    def get_aromatic_rings(self):
        return [ring for ring in self.mol.GetRingInfo().AtomRings()
                if self.mol.GetAtomWithIdx(ring[0]).GetIsAromatic()]

    def get_centroid(self, atom_idx):
        atoms = [self.mol.GetAtomWithIdx(a) for a in atom_idx]
        return centroid_coords(atoms)

    def get_normal(self, ring):
        centroid = self.get_centroid(ring)
        coords1 = coords(self.mol.GetAtomWithIdx(ring[0])) - centroid
        coords2 = coords(self.mol.GetAtomWithIdx(ring[1])) - centroid

        normal = np.cross(coords1, coords2)
        normal /= np.linalg.norm(normal)
        return normal

    def init_hbond(self):
        donors = [atom for atom in self.mol.GetAtoms() if self._is_donor(atom)]
        acceptors = [atom for atom in self.mol.GetAtoms() if self._is_acceptor(atom)]
        return donors, acceptors

    def _is_donor(self, atom):
        if atom.GetAtomicNum() in [7, 8]:
            if _get_bonded_hydrogens(atom):
                return True
        return False

    def _is_acceptor(self, atom):
        if atom.GetAtomicNum() == 8:
            return True
        if atom.GetAtomicNum() == 7 and atom.GetExplicitValence() < 4:
            return True
        return False

    def init_saltbridge(self):
        charged = [atom for atom in self.mol.GetAtoms()
                   if atom.GetFormalCharge() != 0]
        if self.is_protein:
            charge_groups = self._symmetric_charged_protein_atoms()
        else:
            charge_groups = self._symmetric_charged_ligand_atoms()
        return charged, charge_groups

    def _symmetric_charged_protein_atoms(self):
        protein_groups = {}
        for protein_atom in self.mol.GetAtoms():
            if atomname(protein_atom) in ['OD1', 'OD2', 'OE1', 'OE2', 'NH1', 'NH2']:
                if resname(protein_atom) not in protein_groups:
                    protein_groups[resname(protein_atom)] = []
                protein_groups[(resname(protein_atom))] += [protein_atom]
        return protein_groups

    def _symmetric_charged_ligand_atoms(self):
        ligand_groups = {}
        smartss = [('[CX3](=O)[O-]', 2, [1, 2]),
                   ('[CX3](=[NH2X3+])[NH2X3]', 1, [1, 2])]

        idx_to_atom = {atom.GetIdx(): atom for atom in self.mol.GetAtoms()}

        for smarts, k, v in smartss:
            mol = MolFromSmarts(smarts)
            matches = self.mol.GetSubstructMatches(mol)
            for match in matches:
                ligand_groups[match[k]] = [idx_to_atom[match[_v]] for _v in v]
        return ligand_groups

################################################################################
# Compute atom-level interactions

def _get_bonded_hydrogens(atom):
    hydrogens = []
    for bond in atom.GetBonds():
        if bond.GetBeginAtomIdx() != atom.GetIdx():
            hydrogen = bond.GetBeginAtom()
        else:
            hydrogen = bond.GetEndAtom()
            
        if hydrogen.GetAtomicNum() == 1:
            hydrogens += [hydrogen]
    return hydrogens

def _hbond_hydrogen_angle(acceptor, donor):
    best_angle, best_hydrogen = 0, None
    for hydrogen in _get_bonded_hydrogens(donor):
        _angle = angle_atom(donor, hydrogen, acceptor)
        if _angle > best_angle:
            best_angle = _angle
            best_hydrogen = hydrogen
    return best_hydrogen, best_angle

def _hbond_compute(donor_mol, acceptor_mol, settings, protein_is_donor):
    hbonds = []
    for donor in donor_mol.hbond_donors:
        for acceptor in acceptor_mol.hbond_acceptors:
            for hydrogen in _get_bonded_hydrogens(donor):
                dist = distance(acceptor, hydrogen)
                if dist > settings['hbond_dist_cut']: continue
                angle = angle_atom(donor, hydrogen, acceptor)
                if angle < settings['hbond_angle_cut']: continue

                if protein_is_donor:
                    label = 'hbond_donor'
                    protein_atom = donor
                    ligand_atom = acceptor
                else:
                    label = 'hbond_acceptor'
                    protein_atom = acceptor
                    ligand_atom = donor

                hbonds += [{'label': label,
                            'protein_res': resname(protein_atom),
                            'protein_atom': atomname(protein_atom),
                            'ligand_atom': atomname(ligand_atom),
                            'dist': dist,
                            'angle': angle,
                            'hydrogen': atomname(hydrogen)}]
    return hbonds

def hbond_compute(protein, ligand, settings):
    donor = _hbond_compute(protein, ligand, settings, True)
    acceptor = _hbond_compute(ligand, protein, settings, False)
    return acceptor + donor

def saltbridge_compute(protein, ligand, settings):
    # Note that much of the complexity here stems from taking into account
    # symetric atoms. Specifically for carboxylate and guanidinium groups,
    # we consider not just the atom that is arbitrarily assigned a formal
    # charge, but also the atom that is charged in the other resonance
    # structure.

    saltbridges = []
    for protein_atom in protein.charged:
        for ligand_atom in ligand.charged:
            lig_charge = ligand_atom.GetFormalCharge()
            protein_charge = protein_atom.GetFormalCharge()
            if lig_charge * protein_charge >= 0: continue

            # Expand protein_atom and ligand_atom to all symetric atoms
            # ... think carboxylates and guanidiniums.
            if ('saltbridge_resonance' in settings and
                ligand_atom.GetIdx() in ligand.charge_groups):
                ligand_atoms = ligand.charge_groups[ligand_atom.GetIdx()]
            else:
                ligand_atoms = [ligand_atom]
            
            if ('saltbridge_resonance' in settings and
                resname(protein_atom) in protein.charge_groups):
                protein_atoms = protein.charge_groups[resname(protein_atom)]
            else:
                protein_atoms = [protein_atom]

            # Get minimum distance between any pair of protein and ligand
            # atoms in the groups.
            dist = float('inf')
            for _ligand_atom in ligand_atoms:
                for _protein_atom in protein_atoms:
                    _dist = distance(_protein_atom, _ligand_atom)
                    if _dist < dist:
                        dist = _dist
                        closest_protein_atom = _protein_atom
                        closest_ligand_atom = _ligand_atom

            if dist < settings['sb_dist_cut']:
                saltbridges += [{'label': 'saltbridge',
                                 'protein_res': resname(closest_protein_atom),
                                 'protein_atom': atomname(closest_protein_atom),
                                 'ligand_atom': atomname(closest_ligand_atom),
                                 'dist': dist}]
    return saltbridges

def contact_compute(protein, ligand, settings):
    protein = protein.contacts
    ligand = ligand.contacts

    dists = protein[0].reshape(1, -1, 3) - ligand[0].reshape(-1, 1, 3)
    dists = np.linalg.norm(dists, axis=2)
    vdw = protein[1].reshape(1, -1) + ligand[1].reshape(-1, 1)
    contact_idx = np.argwhere(dists < vdw*settings['contact_scale_cut'])

    contacts = []
    for i, j in contact_idx:
        contacts += [{'label': 'contact',
                      'protein_res': protein[2][j],
                      'protein_atom': protein[3][j],
                      'ligand_atom': ligand[3][i],
                      'dist': dists[i, j],
                      'vdw': vdw[i, j]}]
    return contacts

def pipi_compute(protein, ligand, settings):
    pipis = []
    for prot_centroid, prot_normal, prot_res, prot_atom in zip(*protein.pipi):
        for lig_centroid, lig_normal, lig_res, lig_atom in zip(*ligand.pipi):
            displacement = prot_centroid-lig_centroid
            dist = np.linalg.norm(displacement)

            n1_n2 = angle_vector(prot_normal, lig_normal)
            n1_centroid = angle_vector(prot_normal, displacement)
            n2_centroid = angle_vector(lig_normal, displacement)

            # Pi stack
            if (dist < settings['pipi_dist_cut']
                and n1_n2 < settings['pipi_norm_norm_angle_cut'] 
                and n1_centroid < settings['pipi_norm_centroid_angle_cut']
                and n2_centroid < settings['pipi_norm_centroid_angle_cut']):
                pipis += [{'label': 'pipi',
                           'protein_res': prot_res,
                           'protein_atom': prot_atom,
                           'ligand_atom': lig_atom,
                           'dist': dist}]

            # T stack
            elif (dist < settings['pipi_t_dist_cut']
                and settings['pipi_t_norm_norm_angle_cut'] < n1_n2
                and (min(n1_centroid, n2_centroid) < settings['pipi_t_norm_centroid_angle_cut'])):
                pipis += [{'label': 'pi-t',
                           'protein_res': prot_res,
                           'protein_atom': prot_atom,
                           'ligand_atom': lig_atom,
                           'dist': dist}]
    return pipis

################################################################################
# Compute residue-level scores.

def _piecewise(data, opt, cut):
    slope = 1 / (cut-opt)
    intercept = cut * slope

    data = intercept - slope * data
    data[data > 1] = 1
    data[data < 0] = 0
    return data

def _groupby_subset(df, index, col):
    return df[index+[col]].groupby(index)

def nodigits(s):
    return ''.join([i for i in s if not i.isdigit()])

def compute_scores(raw, settings):
    if settings['level'] == 'atom':
        raw['protein_res'] = [r['protein_res']+':'+nodigits(r['protein_atom'])
                              for _, r in raw.iterrows()]

    scores = []
    for label, group in raw.groupby('label'):
        group = group.copy()
        if label == 'pipi':
            group['score'] = _piecewise(group['dist'],
                                        settings['pipi_dist_opt'],
                                        settings['pipi_dist_cut'])
        elif label == 'pi-t':
            group['score'] = _piecewise(group['dist'],
                                        settings['pipi_t_dist_opt'],
                                        settings['pipi_t_dist_cut'])
        elif label == 'contact':
            group['score'] = _piecewise(group['dist'] / group['vdw'],
                                        settings['contact_scale_opt'],
                                        settings['contact_scale_cut'])
        elif label == 'saltbridge':
            group['score'] = _piecewise(group['dist'],
                                        settings['sb_dist_opt'],
                                        settings['sb_dist_cut'])
        elif label in ['hbond_donor', 'hbond_acceptor']:
            group['score'] = (  _piecewise(group['dist'],
                                           settings['hbond_dist_opt'],
                                           settings['hbond_dist_cut'])
                              * _piecewise(180 - group['angle'],
                                           settings['hbond_angle_opt'],
                                           settings['hbond_angle_cut']))

            # One hydrogen bond per hydrogen
            if label == 'hbond_donor':
                idx = _groupby_subset(group,
                                      ['pose', 'protein_res', 'hydrogen'],
                                      'score').idxmax()
            else:
                idx = _groupby_subset(group,
                                      ['pose', 'hydrogen'],
                                      'score').idxmax()
            idx = idx['score']
            group = group.loc[idx]
        group = _groupby_subset(group, ['pose', 'label', 'protein_res'], 'score').sum()
        scores += [group]
    return pd.concat(scores).sort_index()

################################################################################

def fingerprint(protein, ligand, settings):
    fp  = hbond_compute(protein, ligand, settings)
    fp += saltbridge_compute(protein, ligand, settings)
    fp += contact_compute(protein, ligand, settings)
    fp += pipi_compute(protein, ligand, settings)
    return pd.DataFrame.from_dict(fp)

def fingerprint_poseviewer(input_file, poses, settings):
    fps = []
    with gzip.open(input_file) as fp:
        mols =  MaeMolSupplier(fp, removeHs=False)
        protein = Molecule(next(mols), True, settings)
        
        for i, ligand in enumerate(mols):
            if i == poses: break
            if ligand is None:
                print('ligand unreadable')
                continue

            ligand = Molecule(ligand, False, settings)
            fps += [fingerprint(protein, ligand, settings)]
            fps[-1]['pose'] = i

    fps = pd.concat(fps, ignore_index=True, sort=False)
    if 'hydrogen' not in fps:
        fps['hydrogen'] = ''
    fps.loc[fps['hydrogen'].isna(), 'hydrogen'] = ''
    return fps

def ifp(settings, input_file, output_file, poses, convert=False):
    settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}

    if convert:
        temp = tempfile.NamedTemporaryFile(suffix='.maegz')
        convert_mae(input_file, temp.name, poses)
        input_file = temp.name

    # Compute atom-level interactions.
    fps = fingerprint_poseviewer(input_file, poses, settings)

    # Compute residue-level scores.
    scores = compute_scores(fps, settings)

    # Write to files
    fps = fps.set_index(['pose', 'label', 'protein_res', 'protein_atom', 'ligand_atom'])
    fps = fps.sort_index()
    base = output_file.split('.')
    base, ext = base[:-1], base[-1]
    raw_file = '.'.join(base) + '_raw.' + ext

    fps.to_csv(raw_file)
    scores.to_csv(output_file)

@click.command()
@click.argument('input_file')
@click.argument('output_file')
@click.argument('poses', default=100)
@click.option('--convert', is_flag=True)
@click.option('--level', default='residue')
@click.option('--hbond_dist_cut', default=3.0)
@click.option('--hbond_dist_opt', default=2.5)
@click.option('--hbond_angle_cut', default=90.0)
@click.option('--hbond_angle_opt', default=60.0)
@click.option('--sb_dist_cut', default=5.0)
@click.option('--sb_dist_opt', default=4.0)
@click.option('--contact_scale_cut', default=1.75)
@click.option('--contact_scale_opt', default=1.50)
@click.option('--pipi_dist_cut', default=8.0)
@click.option('--pipi_dist_opt', default=7.0)
@click.option('--pipi_norm_norm_angle_cut',  default=30.0)
@click.option('--pipi_norm_centroid_angle_cut', default=45.0)
@click.option('--pipi_t_dist_cut', default=6.0)
@click.option('--pipi_t_dist_opt', default=5.0)
@click.option('--pipi_t_norm_norm_angle_cut', default=60.0)
@click.option('--pipi_t_norm_centroid_angle_cut', default=45.5)
def main(input_file, output_file, poses, convert, **settings):
    ifp(settings, input_file, output_file, poses, convert)

if __name__ == '__main__':
    main()
