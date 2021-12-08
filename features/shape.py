import subprocess
import os
import tempfile
import numpy as np
from schrodinger.structure import StructureReader, StructureWriter

CMD = '$SCHRODINGER/shape_screen -shape {poses1} -screen {poses2} {typing} {norm} -distinct -inplace -NOJOBID'

def write_and_name(conformers, fname):
    titles = []
    with StructureWriter(fname) as writer:
        for i, st in enumerate(conformers):
            st = st.copy()
            assert '-conf-' not in st.title
            st.title = st.title + '-conf-{}'.format(i)
            writer.append(st)
            titles += [st.title]
    return titles

def shape(conformers1, conformers2, version='pharm_max'):
    typing, norm = version.split('_')

    if typing == 'pharm':
        typing = '-pharm'
    elif typing == 'mmod':
        typing = '-atomtypes mmod'
    elif typing == 'element':
        typing = '-atomtypes element'
    elif typing == 'qsar':
        typing = '-atomtypes qsar'
    else:
        assert False, 'Typing {} not supported.'.format(typing)

    if norm == 'max':
        norm = '-norm 1'
    elif norm == 'min':
        norm = '-norm 2'
    else:
        assert False, 'Norm {} not supported.'.format(norm)

    with tempfile.TemporaryDirectory() as wd:
        poses1 = wd+'/poses1.maegz'
        poses2 = wd+'/poses2.maegz'
        output = wd+'/poses1_align.maegz'
        log    = wd+'/poses1_shape.log'

        ligands1 = write_and_name(conformers1, poses1)
        ligands2 = write_and_name(conformers2, poses2)


        cmd = CMD.format(poses1=os.path.basename(poses1),
                         poses2=os.path.basename(poses2),
                         typing=typing, norm=norm)
        print(cmd)
        subprocess.run(cmd, shell=True, cwd=wd)
        print('subprocess complete')

        if not os.path.exists(output):
            with open(log) as fp:
                txt = fp.read()
            assert 'Reference shape must contain at least 3 spheres' in txt, txt
            return 0.5*np.ones((len(ligands1), len(ligands2)))

        sims = np.zeros((len(ligands1), len(ligands2)))
        with StructureReader(output) as sts:
            for k, st in enumerate(sts):
                i = k % len(ligands1)
                j = int(st.title.split('-conf-')[-1])
                sims[i, j] = st.property['r_phase_Shape_Sim']
    return sims
