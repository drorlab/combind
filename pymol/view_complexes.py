from pymol import cmd
from glob import glob

def load_complexes(protein, n=1):
    n = int(n)
    for prot in sorted(glob('{}/structures/proteins/*_prot.mae'.format(protein)))[:n]:
        pdb = prot.split('/')[-1].split('_')[0]
        load_crystal_protein(protein, pdb)
        load_crystal_pose(protein, pdb)
        
    cmd.util.cbao("prot_*")
    cmd.util.cbay("het and crystal_*")
    cmd.show('sticks', "het and crystal_*")
    cmd.hide('lines', 'element h')
    cmd.show('spheres', 'het and prot* and (crystal* expand 5)')

    cmd.show('cartoon')
    cmd.set('cartoon_oval_length', '0.5')
    cmd.set('cartoon_transparency', '0.5')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')
    cmd.hide('everything', 'name H')

def load_crystal_protein(protein, ligand):
    cmd.load('{}/structures/proteins/{}_prot.mae'.format(protein, ligand))
    cmd.set_name('{}_prot'.format(ligand), 'prot_{}'.format(ligand))

def load_crystal_pose(protein, ligand):
    cmd.load('{}/structures/ligands/{}_lig.mae'.format(protein, ligand, ligand))
    cmd.set_name('{}_lig'.format(ligand), 'crystal_{}'.format(ligand))

###################################################################
def load_pose(protein, ligand, struct, pose, prefix):
    pv = glob('{}/docking/{}/{}*{}/*pv.maegz'.format(protein, 'confgen_es4', ligand, struct))[0]
    name = pv.split('/')[-1].split('.')[0]
    
    cmd.load(pv)
    cmd.ungroup('*')
    if pose == 0:
        obj = '{}.{}_lig'.format(name, ligand)
    elif pose < 10:
        obj = '{}.{}_lig_0{}'.format(name, ligand, pose)
    else:
        obj = '{}.{}_lig_{}'.format(name, ligand, pose)
    cmd.set_name(obj, '{}_{}'.format(prefix, ligand))
    cmd.delete(name + '*')

def load_top_glide(protein, n = 1):
    n = int(n)
    grid = None
    for prot in sorted(glob('{}/structures/proteins/*_prot.mae'.format(protein)))[:n]:
        pdb = prot.split('/')[-1].split('_')[0]
        if grid is None:
            grid = pdb
        print(pdb)
        load_pose(protein, pdb, grid, 0, 'glide')
    cmd.show('sticks', "glide_*")
    cmd.hide('lines', 'element h')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')

def load_results(protein, scores):
    struct = glob('{}/docking/grids/*'.format(protein))[0].split('/')[-1]
    load_crystal_protein(protein, struct)
    with open(scores) as fp:
        fp.readline()
        for line in fp:
            if line[:3] == 'com': continue
            (ligand,
             combind_rank, combind_rmsd,
             glide_rank, glide_rmsd,
             best_rank, best_rmsd) = line.strip().split(',')
            ligand = ligand.replace('_lig', '')
            load_pose(protein, ligand, struct, int(combind_rank), 'combind')
            load_pose(protein, ligand, struct, int(glide_rank), 'glide')
            if ligand[:6] != 'CHEMBL':
                load_crystal_pose(protein, ligand)

    cmd.show('sticks', "glide_*")
    cmd.show('sticks', "combind_*")
    cmd.show('sticks', "crystal_*")
    cmd.hide('lines', 'element h')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')

    cmd.util.cbaw('*')
    cmd.color('yellow', 'glide* and element c')
    cmd.color('cyan', 'combind* and element c')
    cmd.set('stick_radius', '0.13')
    
    
    cmd.show('cartoon')
    cmd.set('cartoon_oval_length', '0.5')
    cmd.set('cartoon_transparency', '0.5')

###############################################################

def parse_fp_file(fp_file):
    ifps = {}
    try:
        with open(fp_file) as f:
            pose_num = 0
            for line in f:
                if line.strip() == '': continue
                if line[:4] == 'Pose':
                    pose_num = int(line.strip().split(' ')[1])
                    ifps[pose_num] = {}
                    continue
                sc_key, sc = line.strip().split('=')
                i,r,ss = sc_key.split('-')
                i = int(i)
                sc = float(sc)
                prev_sc = ifps[(i, r)] if (i,r) in ifps[pose_num] else 0
                ifps[pose_num][(i,r)] = max(prev_sc, sc)

    except Exception as e:
        print(e)
        print(fp_file, 'fp not found')
    if len(ifps) == 0:
        print('check', fp_file)
        return {}
    return ifps

def show_interactions(protein, ligand, struct, ifp, pose):
    ifp_file = '{}/ifp/{}/{}_lig-to-{}-confgen_es4.fp'.format(protein, ifp, ligand, struct)
    print(ifp_file)
    ifp = parse_fp_file(ifp_file)[int(pose)]
    cmd.hide('labels')
    cmd.set('label_size', 50)
    for (i, r), score in ifp.items():
        if i not in [2, 3]: continue
        if score < 0.5: continue
        res = r.split(':')[1].split('(')[0]
        cmd.label('{}/ca'.format(res), score)
    
cmd.extend('load_complexes', load_complexes)
cmd.extend('load_top_glide', load_top_glide)
cmd.extend('load_results', load_results)
cmd.extend('show_interactions', show_interactions)
