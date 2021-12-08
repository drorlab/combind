from pymol import cmd
import sys

def read_results(fname):
    data = {}
    with open(fname) as fp:
        fp.readline()
        for line in fp:
            if line[:3] == 'com': continue
            print(line.strip().split(','))
            ligand, combind, _, glide, _, best, _ = line.strip().split(',')

            data[ligand] = (int(combind),
                            int(glide),
                            int(best) if best != 'None' else None)
    return data

def pose(poseviewer, pose_number, label = 'pose'):
    name = '.'.join(poseviewer.split('/')[-1].split('.')[:-1])
    struct = name.split('-to-')[-1].split('_')[0] + '_prot'
    ligand = name.split('-to-')[0]
    if pose_number == 0:
        pose_number = ''
    elif pose_number < 10:
        pose_number = '0' + str(pose_number)
    else:
        pose_number = str(pose_number)
    print(name, struct, ligand, pose_number)
    cmd.load(poseviewer)
    cmd.split_states(name)
    cmd.delete(name)

    pose_name = label + '_' + ligand
    grid_name = 'grid_' + ligand
    cmd.set_name(ligand + pose_number, pose_name)
    cmd.set_name(struct, grid_name)
    cmd.delete(ligand + '*')

    cmd.show_as('cartoon', grid_name)
    cmd.show_as('sticks', pose_name)

def results(fname, protein, struct, docking):
    data = read_results(fname)
    for ligand, (combind, glide, best) in data.items():
        print(ligand)
        poseviewer = '{}/docking/{}/{}-to-{}/{}-to-{}_pv.maegz'.format(protein, docking,
                                                                       ligand, struct,
                                                                       ligand, struct)
        pose(poseviewer, combind, 'combind')
        pose(poseviewer, glide,   'glide')

    cmd.util.cbam('combind*')
    cmd.util.cbag('glide*')
    cmd.util.cbaw('grid*')
    cmd.hide('sticks', 'element h and not (element o+n extend 1)')

    cmd.group('glide', 'glide*')
    cmd.group('combind', 'combind*')
    cmd.group('grid', 'grid*')


cmd.extend('pose', pose)
cmd.extend('results', results)
