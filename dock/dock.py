from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.rmsd import ConformerRmsd
import os
import subprocess

GLIDE_ES4 = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   SP
NENHANCED_SAMPLING   4
'''

GLIDE = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   30
POSTDOCK_NPOSE   30
PRECISION   SP
'''

def docking_failed(glide_log):
    if not os.path.exists(glide_log):
        return False
    with open(glide_log) as fp:
        logtxt = fp.read()
    phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
               'NO VALID POSES AFTER MINIMIZATION: SKIPPING.',
               'No Ligand Poses were written to external file',
               'GLIDE WARNING: Skipping refinement, etc. because rough-score step failed.']
    return any(phrase in logtxt for phrase in phrases)

def dock(grid, ligands, root, name, enhanced, infile=None, reference=None):
    if infile is None:
        infile = GLIDE_ES4 if enhanced else GLIDE
    glide_in = '{}/{}.in'.format(root, name)
    glide_pv = '{}/{}_pv.maegz'.format(root, name)
    glide_log = '{}/{}.log'.format(root, name)
    glide_cmd = 'glide -WAIT -LOCAL -RESTART {}'.format(os.path.basename(glide_in))

    if os.path.exists(glide_pv):
        return

    if enhanced and docking_failed(glide_log):
        return

    if not os.path.exists(root):
        os.system('mkdir {}'.format(root))
    with open(glide_in, 'w') as fp:
        fp.write(infile.format(grid=grid, ligands=ligands, reference=reference))

    subprocess.run(glide_cmd, cwd=root, shell=True)

def filter_native(native, pv, out, thresh):
    with StructureReader(native) as sts:
        native = list(sts)
        assert len(native) == 1, len(native)
        native = native[0]

    near_native = []
    with StructureReader(pv) as reader:
        receptor = next(reader)
        for st in reader:
            conf_rmsd = ConformerRmsd(native, st)
            if conf_rmsd.calculate() < thresh:
                near_native += [st]

    print('Found {} near-native poses'.format(len(near_native)))
    if not near_native:
        print('Resorting to native pose.')
        native.property['r_i_docking_score'] = -10.0
        near_native = [native]

    with StructureWriter(out) as writer:
        writer.append(receptor)
        for st in near_native:
            writer.append(st)
