from pymol import cmd
import sys
import pandas as pd

def style():
	cmd.show('cartoon')
	cmd.show('lines')
	cmd.hide('sticks')
	cmd.util.cbaw()
	cmd.color('slate', 'het and element C')
	cmd.hide('everything', 'element H and (element C extend 1)')

def pose_name(group, pose):
	if pose == 0:
		pose = group.split('-')[0]
	elif pose < 10:
		pose = '_0{}'.format(pose)
	else:
		pose = '_{}'.format(pose)
	return '{}*.*{}'.format(group, pose)

def enable(group, pose, prot=True):
	cmd.enable(group+'*-to-*')
	cmd.disable(group+'*-to-*.*')
	cmd.enable(pose_name(group, pose))
	if prot:
		cmd.enable('{}*.*_prot'.format(group))

def show_interactions(ifp_file, interaction, lig, pose, delete=True, disable=True):
	pose = int(pose)

	if delete:
		cmd.delete('dist*')
		cmd.delete('ps*')
	if disable:
		cmd.disable('*')
	style()

	enable(lig, pose)

	if interaction == 'all':
		for interaction in ['sb', 'hbond', 'contact', 'pipi']:
			 show_interactions(ifp_file, interaction, lig, pose,
			                   delete=False, disable=False)
	
	df = pd.read_csv(ifp_file)

	if interaction == 'hbond':
		idx = df['label'] == 'hbond_acceptor'
		idx |= df['label'] == 'hbond_donor'
		thresh = 3.5
		color = 'yellow'
	elif interaction == 'sb':
		idx = df['label'] == 'saltbridge'
		thresh = 4
		color = 'magenta'
	elif interaction == 'contact':
		idx = df['label'] == 'contact'
		thresh = 1.25
		color='smudge'
	elif interaction == 'pipi':
		idx = df['label'] == 'pipi'
		idx |= df['label'] == 'pi-t'
		thresh = 7.0
		color='green'

	idx &= df['pose'] == pose

	for i, row in df[idx].iterrows():
		if interaction == 'contact' and row['dist'] > thresh*row['vdw']: continue
		if interaction != 'contact' and row['dist'] > thresh: continue
		
		chain, resid, _, _ = row['protein_res'].split(':')
		prot = '{}*.*prot and chain {} and resid {} and name {}'.format(lig,
		                                                               chain,
		                                                               resid,
		                                                               row['protein_atom'].replace(',', '+'))
		ligand = '{} and name {}'.format(pose_name(lig, pose), row['ligand_atom'].replace(',', '+'))

		print(prot, ligand)

		cmd.pseudoatom('ps{}{}prot'.format(interaction, i), prot)
		cmd.pseudoatom('ps{}{}lig'.format(interaction, i), ligand)

		cmd.dist('dist'+interaction+str(i),
		         'ps{}{}prot'.format(interaction, i),
		         'ps{}{}lig'.format(interaction, i))
		cmd.color(color, 'dist'+interaction+str(i))

	cmd.set('dash_width', 6)
	cmd.set('dash_width', 3, 'distcontact*')

	cmd.enable('{}*.*prot'.format(lig))
	cmd.enable(pose_name(lig, pose))
	cmd.enable(lig)

cmd.extend('show_interactions', show_interactions)
