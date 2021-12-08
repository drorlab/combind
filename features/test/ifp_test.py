import pytest
import ifp
import gzip
from rdkit.Chem.rdmolfiles import MaeMolSupplier


settings = {'version'           : 'rd1',
            'level'             : 'residue',
            'hbond_dist_opt'    : 2.5,
            'hbond_dist_cut'    : 3.0,
            'hbond_angle_opt'   : 60.0,
            'hbond_angle_cut'   : 90.0,
            'sb_dist_opt'       : 4.0,
            'sb_dist_cut'       : 5.0,
            'contact_scale_opt' : 1.25,
            'contact_scale_cut' : 1.75,
            'pipi_dist_cut'     : 7.0,
            'pipi_dist_opt'     : 7.0,
            'pipi_norm_norm_angle_cut'     : 30.0,
            'pipi_norm_centroid_angle_cut' : 45.0,
            'pipi_t_dist_cut': 6.0,
            'pipi_t_dist_opt': 5.0,
            'pipi_t_norm_norm_angle_cut': 60.0,
            'pipi_t_norm_centroid_angle_cut': 45.0}

settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}

with gzip.open('test/3ZPR_lig-to-2VT4_pv.maegz') as fp:
    mols =  MaeMolSupplier(fp, removeHs=False)
    protein = ifp.Molecule(next(mols), True, settings)
    ligands = [ifp.Molecule(mol, False, settings) for mol in mols]

settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
def test_version():
    import rdkit
    assert rdkit.__version__ == '2020.03.1'

def test_hydrogenbond():
    i = ifp.hbond_compute(protein, ligands[0], settings)
    assert len(i) == 2

def test_saltbridge_none():
    i = ifp.saltbridge_compute(protein, ligands[0], settings)
    assert len(i) == 0

def test_saltbridge_one():
    i = ifp.saltbridge_compute(protein, ligands[3], settings)
    assert len(i) == 1

def test_contact():
    ifp.contact_compute(protein, ligands[0], settings)

def test_pipi_tstack():
    i = ifp.pipi_compute(protein, ligands[0], settings)
    assert len(i) == 4
    i = ifp.pipi_compute(protein, ligands[3], settings)
    assert len(i) == 4

def test_pipi_pstack():
    i = ifp.pipi_compute(protein, ligands[173], settings)
    print(i)
    assert len(i) == 4

    i = ifp.pipi_compute(protein, ligands[180], settings)
    print(i)
    assert len(i) == 2

def test_pipi_tstack_6IBL():

    with gzip.open('test/6IBL-to-2VT4_pv.maegz') as fp:
        mols =  MaeMolSupplier(fp, removeHs=False)
        protein = ifp.Molecule(next(mols), True, settings)
        ligands = [ifp.Molecule(mol, False, settings) for mol in mols]
    i = ifp.pipi_compute(protein, ligands[3], settings)
    assert len(i) == 2

#ifp.fingerprint_poseviewer('test/pv.maegz', 100, settings)
