import tempfile
import numpy as np
import subprocess
import os
from schrodinger.structure import StructureReader, StructureWriter, SmilesStructure
from schrodinger.structutils.rmsd import ConformerRmsd
from schrodinger.structutils.analyze import generate_smiles
from schrodinger.structutils.analyze import evaluate_smarts_canvas

def mcss(sts1, sts2, mcss_types_file):
    """
    Computes rmsd between mcss for atoms in two poseviewer files.

    Returns a (# poses in pv1) x (# poses in pv2) np.array of rmsds.
    """
    memo = {}
    sts1 = [merge_halogens(st.copy()) for st in sts1]
    sts2 = [merge_halogens(st.copy()) for st in sts2]

    rmsds = []
    for j, st1 in enumerate(sts1):
        smi1 = generate_smiles(st1)
        n_st1_atoms = n_atoms(st1)
        rmsds += [np.zeros(len(sts2)) + float('inf')]
        for i, st2 in enumerate(sts2):
            smi2 = generate_smiles(st2)
            n_st2_atoms = n_atoms(st2)

            if (smi1, smi2) in memo:
                mcss, n_mcss_atoms = memo[(smi1, smi2)]
            else:
                mcss = compute_mcss(st1, st2, mcss_types_file)
                # Capitalizing is useful to prevent invalid structures.
                # This structure isn't used for anything other than
                # counting atoms.
                mcss_st = SmilesStructure(mcss['st1'][0].split(',')[0].upper()).get2dStructure()
                n_mcss_atoms = n_atoms(mcss_st)
                memo[(smi1, smi2)] = (mcss, n_mcss_atoms)
                memo[(smi2, smi1)] = ({'st1': mcss['st2'], 'st2': mcss['st1']}, n_mcss_atoms)

            if (2*n_mcss_atoms <= min(n_st1_atoms, n_st2_atoms)
                or n_mcss_atoms <= 10):
                continue

            rmsds[-1][i] = compute_mcss_rmsd(st1, st2, mcss)

    return np.vstack(rmsds)

def compute_mcss_rmsd(st1, st2, mcss):
    """
    Compute minimum rmsd between mcss(s).

    Takes into account that there can be multiple mcss smarts patterns (
    i.e. two patterns that are the same size) and each smarts pattern could
    map to multiple atom indices (e.g. symetric groups).
    """
    rmsd = float('inf')
    for smarts1, smarts2 in zip(mcss['st1'], mcss['st2']):
        atom_idxs1 = evaluate_smarts_canvas(st1, smarts1)
        atom_idxs2 = evaluate_smarts_canvas(st2, smarts2)

        for atom_idx1 in atom_idxs1:
            for atom_idx2 in atom_idxs2:
                _rmsd = calculate_rmsd(st1, st2, atom_idx1, atom_idx2)
                rmsd = min(_rmsd, rmsd)
    return rmsd

def compute_mcss(st1, st2, mcss_types_file):
    """
    Compute smarts patterns for mcss(s) between two structures.
    """
    cmd = "$SCHRODINGER/utilities/canvasMCS -imae {} -ocsv {} -stop 10 -atomtype C {}"
    with tempfile.TemporaryDirectory() as wd:
        mae = wd+'/temp.maegz'
        csv = wd+'/temp.csv'

        st1.title = 'st1'
        st2.title = 'st2'
        stwr = StructureWriter(mae)
        stwr.append(st1)
        stwr.append(st2)
        stwr.close()

        r = subprocess.run(cmd.format(os.path.basename(mae), os.path.basename(csv),
                                      os.path.abspath(mcss_types_file)),
                           cwd=wd, shell=True, stderr=subprocess.PIPE)


        # mcss can fail with memory usage error, generally when macrocycles
        # are present. just skip such cases.
        failed = 'memory usage' in str(r.stderr)

        if not failed:
            assert os.path.exists(csv)
            mcss = {'st1': [], 'st2': []}
            with open(csv) as fp:
                fp.readline()
                for line in fp:
                    lig = line.strip().split(',')[1]
                    smarts = line.strip().split(',')[-1]
                    mcss[lig] += [smarts]
        else:
            print('mcss failed')
            mcss = {'st1': ['C'], 'st2': ['C']}
    return mcss

def calculate_rmsd(pose1, pose2, atom_idx1, atom_idx2):
    """
    Calculates the RMSD between the atoms atom_idx1 in pose1
    and the atoms atom_idx2 in pose2.

    pose1, pose2: schrodinger.structure
    atom_idx1, atom_idx2: [int, ...]
    merge_halogens: If true then change the atomic number of all halogens
                    to 9 (the atomic number of flourine) before computing
                    rmsds. This allows for MCSS that treat all halogens
                    the same.
    """
    substructure1 = pose1.extract(atom_idx1)
    substructure2 = pose2.extract(atom_idx2)
    try:
        calc = ConformerRmsd(substructure1, substructure2)
        calc.use_heavy_atom_graph = True
        calc.renumber_structures = True
        rmsd = calc.calculate()
    except:
        # This is necessary because there is a bug in the
        # Schrodinger software that results in incorrect
        # atom indices being used when the heavy_atom_graph
        # is used. That being said, the above is more reliable
        # than the below, so should be tried first.
        calc = ConformerRmsd(substructure1, substructure2)
        calc.use_heavy_atom_graph = True
        calc.renumber_structures = False
        rmsd = calc.calculate()
    return rmsd

def merge_halogens(structure):
    """
    Sets atomic number for all halogens to be that for flourine.
    This enable use of ConformerRmsd for atom typing schemes that
    merge halogens.
    """
    for atom in structure.atom:
        if atom.atomic_number in [9, 17, 35, 53]:
            atom.atomic_number = 9
    return structure

def n_atoms(st):
    return sum(atom.element != 'H' for atom in st.atom)
