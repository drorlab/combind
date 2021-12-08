import pytest
from score.pairs import PosePair
from containers import Pose

features = {'sb': [1], 'hbond': [2], 'contact': [11]}

def pose(rmsd=0.0, gscore=0.0, fp={}):
	return Pose(rmsd, gscore, fp)

def test_correct_both():
	pose1 = pose(rmsd=1.0)
	pose2 = pose(rmsd=1.4)
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)
	
	assert pp.correct() == 1.0

def test_correct_one_1():
	pose1 = pose(rmsd=1.0)
	pose2 = pose(rmsd=2.0+1.0)
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)
	
	assert pp.correct() == 0.0

def test_correct_one_2():
	pose1 = pose(rmsd=2.0+1.0)
	pose2 = pose(rmsd=1.4)
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)
	
	assert pp.correct() == 0.0

def test_correct_one_neither():
	pose1 = pose(rmsd=2.0+1.0)
	pose2 = pose(rmsd=2.0+0.1)
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)
	
	assert pp.correct() == 0.0

def test_get_feature_empty():
	pose1 = pose(fp={})
	pose2 = pose(fp={})
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)

	assert pp.overlap('sb') == 0.0
	assert pp.overlap('hbond') == 0.0
	assert pp.overlap('contact') == 0.0
	assert pp.mcss_score == 0.0

def test_get_feature_single():
	pose1 = pose(fp={(1, 23): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 4.0, features, 2.0)

	assert pp.overlap('sb') == 1.0
	assert pp.overlap('hbond') == 0.0
	assert pp.overlap('contact') == 0.0
	assert pp.mcss_score == 4.0

def test_get_feature_mismatch():
	pose1 = pose(fp={(1, 2): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 1.0, features, 2.0)

	assert pp.overlap('sb') == 0.0
	assert pp.overlap('hbond') == 0.0
	assert pp.overlap('contact') == 0.0
	assert pp.mcss_score == 1.0

def test_get_feature_multiple_of_same_type():
	pose1 = pose(fp={(1, 2): 1.0, (1, 23): 1.0})
	pose2 = pose(fp={(1, 2): 0.0, (1, 23): 1.0})
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)

	assert pp.overlap('sb') == 1.0
	assert pp.overlap('hbond') == 0.0
	assert pp.overlap('contact') == 0.0
	assert pp.mcss_score == 0.0

def test_tanimoto_empty():
	pose1 = pose(fp={})
	pose2 = pose(fp={})
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)

	assert pp.tanimoto('sb') == 1/2
	assert pp.tanimoto('hbond') == 1/2
	assert pp.tanimoto('contact') == 1/2
	assert pp.mcss_score == 0.0

def test_tanimoto_single():
	pose1 = pose(fp={(1, 23): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 4.0, features, 2.0)

	assert pp.tanimoto('sb') == 2/3
	assert pp.tanimoto('hbond') == 1/2
	assert pp.tanimoto('contact') == 1/2
	assert pp.mcss_score == 4.0

def test_tanimoto_mismatch():
	pose1 = pose(fp={(1, 2): 1.0})
	pose2 = pose(fp={(1, 23): 1.0})
	pp = PosePair(pose1, pose2, 1.0, features, 2.0)

	assert pp.tanimoto('sb') == 1/4
	assert pp.tanimoto('hbond') == 1/2
	assert pp.tanimoto('contact') == 1/2
	assert pp.mcss_score == 1.0

def test_tanimoto_multiple_of_same_type():
	pose1 = pose(fp={(1, 2): 1.0, (1, 23): 1.0})
	pose2 = pose(fp={(1, 2): 0.0, (1, 23): 1.0})
	pp = PosePair(pose1, pose2, 0.0, features, 2.0)

	assert pp.tanimoto('sb') == 2/4
	assert pp.tanimoto('hbond') == 1/2
	assert pp.tanimoto('contact') == 1/2
	assert pp.mcss_score == 0.0
