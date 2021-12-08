"""
Tests for the LigPair class.
"""

import pytest
from score.pairs import LigPair
from containers import Ligand, Pose


def create_ligand(name, fps):
	ligand = Ligand(name, {}, {})
	ligand.poses = []
	for fp in fps:
		ligand.poses += [Pose(0, 0, fp)]
	return ligand

def test_empty_tanimoto():
	lig1 = create_ligand('lig1', [{}]*5)
	lig2 = create_ligand('lig2', [{}]*5)

	lp = LigPair(lig1, lig2, {'sb': [1], 'hbond': [2]}, None, 4)
	lp.init_pose_pairs()

	assert len(lp.pose_pairs) == 16

	for i in range(4):
		for j in range(4):
			assert lp.get_feature('sb', i, j) == 1/2
			assert lp.get_feature('hbond', i, j) == 1/2

def test_max_one_tanimoto():
	lig1 = create_ligand('lig1', [{(1, 23): 1.0}, {(1, 23): 0.5}, {}])
	lig2 = create_ligand('lig2', [{}, {(1, 23): 0.5}, {(1, 23): 1.0}])

	lp = LigPair(lig1, lig2, {'sb': [1], 'hbond': [2]}, None, 4)
	lp.init_pose_pairs()

	assert len(lp.pose_pairs) == 9

	assert lp.get_feature('sb', 0, 0) == 1/3
	assert lp.get_feature('sb', 1, 0) == 1/2.5
	assert lp.get_feature('sb', 2, 0) == 1/2
	assert lp.get_feature('sb', 0, 1) == (1+0.5**0.5) / (3.5 - 0.5**0.5)
	assert lp.get_feature('sb', 1, 1) == 1.5/2.5
	assert lp.get_feature('sb', 2, 1) == 1/2.5
	assert lp.get_feature('sb', 0, 2) == 2/3
	assert lp.get_feature('sb', 1, 2) == (1+0.5**0.5) / (3.5 - 0.5**0.5)
	assert lp.get_feature('sb', 2, 2) == 1/3

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) == 1/2

def test_max_less_than_one_tanimoto():
	lig1 = create_ligand('lig1', [{(1, 23): 0.9}, {(1, 23): 0.5}, {}])
	lig2 = create_ligand('lig2', [{}, {(1, 23): 0.5}, {(1, 23): 1.0}])

	lp = LigPair(lig1, lig2, {'sb': [1], 'hbond': [2]}, None, 4)
	lp.init_pose_pairs()

	assert len(lp.pose_pairs) == 9

	assert lp.get_feature('sb', 0, 0) == 1/2.9
	assert lp.get_feature('sb', 1, 0) == 1/2.5
	assert lp.get_feature('sb', 2, 0) == 1/2
	assert lp.get_feature('sb', 0, 1) == (1 + (0.9*0.5)**0.5) / (3.4 - (0.9*0.5)**0.5)
	assert lp.get_feature('sb', 1, 1) == 1.5 / 2.5
	assert lp.get_feature('sb', 2, 1) == 1 / 2.5
	assert lp.get_feature('sb', 0, 2) == (1 + .9**0.5) / (3.9 - .9**0.5)
	assert lp.get_feature('sb', 1, 2) == (1+0.5**0.5) / (3.5 - 0.5**0.5)
	assert lp.get_feature('sb', 2, 2) == 1/3

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) == 1/2

def test_max_greater_than_one_tanimoto():
	lig1 = create_ligand('lig1', [{(1, 23): 1.0, (1, 20): 1.0},
	                     		  {(1, 23): 0.5},
	                     		  {(1, 20): 1.0}])
	lig2 = create_ligand('lig2', [{},
	                     		  {(1, 23): 0.5},
	                     		  {(1, 23): 1.0, (1, 20): 1.0}])

	lp = LigPair(lig1, lig2, {'sb': [1], 'hbond': [2]}, None, 4)
	lp.init_pose_pairs()

	assert len(lp.pose_pairs) == 9

	assert lp.get_feature('sb', 0, 0) == 1/4
	assert lp.get_feature('sb', 1, 0) == 1/2.5
	assert lp.get_feature('sb', 2, 0) == 1/3
	assert lp.get_feature('sb', 0, 1) == (1 + 0.5**0.5) / (4.5 - 0.5**0.5)
	assert lp.get_feature('sb', 1, 1) == 1.5 / 2.5
	assert lp.get_feature('sb', 2, 1) == 1/3.5
	assert lp.get_feature('sb', 0, 2) == 3/4
	assert lp.get_feature('sb', 1, 2) == (1 + 0.5**0.5)/ (4.5 - 0.5**0.5)
	assert lp.get_feature('sb', 2, 2) == 2/4

	for i in range(3):
		for j in range(3):
			assert lp.get_feature('hbond', i, j) == 1/2
