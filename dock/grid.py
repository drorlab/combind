import os
import subprocess
from glob import glob
from schrodinger.structure import StructureReader
from schrodinger.structutils.transform import get_centroid

GRID_IN = """
GRID_CENTER {x},{y},{z}
GRIDFILE {pdb}.zip
INNERBOX 15,15,15
OUTERBOX 30,30,30
RECEP_FILE {prot}
"""

CMD = "glide -WAIT {infile}"
INFILE = '{pdb}.in'
ZIPFILE = '{pdb}.zip'

def centroid(ligfile):
    with StructureReader(ligfile) as st:
        st = next(st)
    c = get_centroid(st)
    x,y,z = c[:3]
    return x, y, z

def make_grid(pdb,
              PROTFILE='structures/proteins/{pdb}_prot.mae',
              LIGFILE='structures/ligands/{pdb}_lig.mae',
              CWD='structures/grids/{pdb}',
              grid_in=None):
    if grid_in is None:
        grid_in = GRID_IN

    cwd = os.path.abspath(CWD.format(pdb=pdb))
    zipfile = os.path.abspath(cwd+'/'+ZIPFILE.format(pdb=pdb))
    infile = os.path.abspath(cwd+'/'+INFILE.format(pdb=pdb))
    ligfile = os.path.abspath(LIGFILE.format(pdb=pdb))
    protfile = os.path.abspath(PROTFILE.format(pdb=pdb))
    cmd = CMD.format(infile=os.path.basename(infile))

    if os.path.exists(zipfile):
        return # Done.
    if not (os.path.exists(ligfile) and os.path.exists(protfile)):
        print(ligfile, protfile)
        return # Not ready.

    print('making grid', pdb)

    for path in glob(cwd + '/*'):
        os.remove(path)
    os.makedirs(cwd, exist_ok=True)

    x, y, z = centroid(ligfile)

    with open(infile, 'w') as fp:
        fp.write(grid_in.format(x=x, y=y, z=z, pdb=pdb, prot=protfile))

    subprocess.run(cmd, cwd=cwd, shell=True)
