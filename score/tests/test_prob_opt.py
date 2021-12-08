"""
Tests for optimization code.
"""

import pytest
import numpy as np

from score.prob_opt import PredictStructs
from score.density_estimate import DensityEstimate
from containers import Ligand, Pose

###############################################################################
# Create PredictStructs objects

def create_ligand(name, rmsds, gscores, fps):
	ligand = Ligand(name, {}, {})
	ligand.poses = []
	for r, g, f in zip(rmsds, gscores, fps):
		ligand.poses += [Pose(r, g, f)]
	return ligand

def basic_ps():
	ligands = {'lig1': create_ligand('lig1', [0, 1], [-2, -1.5], [{}, {(1, 23): 1.0}]),
			   'lig2': create_ligand('lig2', [0, 1], [-10, -3.5], [{(1, 23): 1.0}, {}])}
	stats = {'native': {'sb': DensityEstimate(domain = (0, 1)),
						'hbond': DensityEstimate(domain = (0, 1))},
			 'reference': {'sb': DensityEstimate(domain = (0, 1)),
			 			   'hbond': DensityEstimate(domain = (0, 1))}}
	stats['native']['sb'].x = np.linspace(0, 1, 100)
	stats['reference']['sb'].x = np.linspace(0, 1, 100)
	stats['native']['sb'].fx = np.linspace(0, 2, 100)
	stats['reference']['sb'].fx = np.linspace(2, 0, 100)
	features = {'sb': [1]}
	return PredictStructs(ligands, None, stats, features, 3, 1.0)

def two_interaction_ps():
	ligands = {'lig1': create_ligand('lig1', [0, 1], [-2, -1.5],
	                                 [{(2, 10): 0.0}, {(1, 23): 1.0, (2, 10): 1.0}]),
			   'lig2': create_ligand('lig2', [0, 1], [-10, -3.5],
			                         [{(1, 23): 1.0}, {(2, 10): 0.6}])}
	stats = {'native': {'sb': DensityEstimate(domain = (0, 1)),
						'hbond': DensityEstimate(domain = (0, 1))},
			 'reference': {'sb': DensityEstimate(domain = (0, 1)),
			 			   'hbond': DensityEstimate(domain = (0, 1))}}
	stats['native']['sb'].x = np.linspace(0, 1, 100)
	stats['reference']['sb'].x = np.linspace(0, 1, 100)
	stats['native']['sb'].fx = np.linspace(0.0, 2, 100)
	stats['reference']['sb'].fx = np.linspace(2, 0.0, 100)

	stats['native']['hbond'].x = np.linspace(0, 1, 100)
	stats['reference']['hbond'].x = np.linspace(0, 1, 100)
	stats['native']['hbond'].fx = np.linspace(0, 2, 100)
	stats['reference']['hbond'].fx = np.linspace(1.0, 1.0, 100)
	features = {'sb': [1], 'hbond': [2, 3]}
	return PredictStructs(ligands, None, stats, features, 3, 1.0)

#############################################################################

def test_get():
	ps = basic_ps()
	assert ps._get_physics_score('lig1', 0) == -2.0
	assert ps._get_physics_score('lig1', 1) == -1.5
	assert ps._get_physics_score('lig2', 0) == -10.0
	assert ps._get_physics_score('lig2', 1) == -3.5
	assert ps._num_poses('lig1') == 2
	assert ps._num_poses('lig2') == 2
	with pytest.raises(KeyError):
		ps._num_poses('lig3')

	assert ps._get_feature('sb', 'lig1', 'lig2', 0, 0) == 1/3
	assert ps._get_feature('sb', 'lig1', 'lig2', 0, 1) == 1/2
	assert ps._get_feature('sb', 'lig1', 'lig2', 1, 0) == 2/3
	assert ps._get_feature('sb', 'lig1', 'lig2', 1, 1) == 1/3
	assert ps._get_feature('sb', 'lig2', 'lig1', 0, 0) == 1/3
	assert ps._get_feature('sb', 'lig2', 'lig1', 0, 1) == 2/3
	assert ps._get_feature('sb', 'lig2', 'lig1', 1, 0) == 1/2
	assert ps._get_feature('sb', 'lig2', 'lig1', 1, 1) == 1/3

def test_like():
	ps = basic_ps()
	def like1(feature, p1, p2):
		return ps._likelihoods_for_pair_and_single_feature(feature,
		            {'lig1':p1, 'lig2':p2},'lig1', 'lig2')
	assert like1('sb', 0, 0) == pytest.approx((1/3, 2/3, 2 - 2/3), 0.001)
	assert like1('sb', 0, 1) == pytest.approx((1/2, 1, 1), 0.001)
	assert like1('sb', 1, 0) == pytest.approx((2/3, 2 - 2/3, 2/3), 0.001)
	assert like1('sb', 1, 1) == pytest.approx((1/3, 2/3, 2 - 2/3), 0.001)
	
	# These four should give identical results to above 4.
	#   just checking that reversing ligand names isn't a problem.
	def like2(feature, p1, p2):
		return ps._likelihoods_for_pair_and_single_feature(feature,
		            {'lig1':p1, 'lig2':p2}, 'lig2', 'lig1')
	assert like2('sb', 0, 0) == pytest.approx((1/3, 2/3, 2 - 2/3), 0.001)
	assert like2('sb', 0, 1) == pytest.approx((1/2, 1, 1), 0.001)
	assert like2('sb', 1, 0) == pytest.approx((2/3, 2 - 2/3, 2/3), 0.001)
	assert like2('sb', 1, 1) == pytest.approx((1/3, 2/3, 2 - 2/3), 0.001)

def test_like_none():
	ps = basic_ps()
	def like1(feature, p1, p2):
		return ps._likelihoods_for_pair_and_single_feature(feature,
		            {'lig1':p1, 'lig2':p2},'lig1', 'lig2')
	def like2(feature, p1, p2):
		return ps._likelihoods_for_pair_and_single_feature(feature,
		            {'lig1':p1, 'lig2':p2},'lig2', 'lig1')

	for i in range(2):
		for j in range(2):
			assert like1('hbond', i, j) == (0.0, 1.0, 1.0)
			assert like2('hbond', i, j) == (0.0, 1.0, 1.0)

def test_ratio():
	ps = basic_ps()
	def like(p1, p2):
		return ps._log_likelihood_ratio_pair({'lig1':p1, 'lig2':p2}, 'lig1', 'lig2')
	assert like(0, 0) == pytest.approx(np.log(2/3) - np.log(2 - 2/3), 0.001)
	assert like(0, 1) == pytest.approx(np.log(1) - np.log(1), 0.001)
	assert like(1, 0) == pytest.approx(np.log(2 - 2/3) - np.log(2/3), 0.001)
	assert like(1, 1) == pytest.approx(np.log(2/3) - np.log(2 - 2/3), 0.001)


	def like(p1, p2):
		return ps._log_likelihood_ratio_pair({'lig1':p1, 'lig2':p2}, 'lig2', 'lig1')
	assert like(0, 0) == pytest.approx(np.log(2/3) - np.log(2 - 2/3), 0.001)
	assert like(0, 1) == pytest.approx(np.log(1) - np.log(1), 0.001)
	assert like(1, 0) == pytest.approx(np.log(2 - 2/3) - np.log(2/3), 0.001)
	assert like(1, 1) == pytest.approx(np.log(2/3) - np.log(2 - 2/3), 0.001)

def test_ratio_two():
	ps = two_interaction_ps()

	x = ps._log_likelihood_ratio_pair({'lig1':1, 'lig2':1},'lig1', 'lig2')
	y = ps._log_likelihood_ratio_pair({'lig1':1, 'lig2':1},'lig2', 'lig1')

	assert x == y
	assert x == pytest.approx(np.log(2/3) - np.log(2 - 2/3)
	             + np.log(2*(1 +.6**0.5) / (2 + .6 + 1 - .6**0.5))
	             - np.log(1.0), 0.001)


	x, p_x_n, p_x = ps._likelihoods_for_pair_and_single_feature('hbond',
		            {'lig1':0, 'lig2':1}, 'lig2', 'lig1')

	assert x == 1 / 2.6
	assert p_x_n == 2 * 1 / (2 + .6)
	assert p_x == 1.0

	x = ps._log_likelihood_ratio_pair({'lig1':0, 'lig2':1},'lig1', 'lig2')
	y = ps._log_likelihood_ratio_pair({'lig1':0, 'lig2':1},'lig2', 'lig1')

	assert x == y
	assert x == pytest.approx(np.log(1) - np.log(1)
	                          + np.log(2 * 1 / (2 + .6)) - np.log(1.0), 0.001)

def test_partial():
	ps = basic_ps()
	def like(p1, p2, lig):
		return ps._partial_log_posterior({'lig1':p1, 'lig2':p2}, lig)

	assert like(0, 0, 'lig1') == np.log(2/3) - np.log(2 - 2/3) + 2*1.0
	assert like(0, 1, 'lig1') == np.log(1) - np.log(1) + 2*1.0
	assert like(1, 0, 'lig1') == np.log(2 - 2/3) - np.log(2/3) + 1.5*1.0
	assert like(1, 1, 'lig1') == np.log(2/3) - np.log(2 - 2/3) + 1.5*1.0

	assert like(0, 0, 'lig2') == np.log(2/3) - np.log(2 - 2/3) + 10*1.0
	assert like(0, 1, 'lig2') == np.log(1) - np.log(1) + 3.5*1.0
	assert like(1, 0, 'lig2') == np.log(2 - 2/3) - np.log(2/3) + 10.0*1.0
	assert like(1, 1, 'lig2') == np.log(2/3) - np.log(2 - 2/3) + 3.5*1.0

def test_posterior():
	ps = basic_ps()
	def like(p1, p2):
		return ps.log_posterior({'lig1':p1, 'lig2':p2})

	assert like(0, 0) == np.log(2/3) - np.log(2 - 2/3) + 2*1.0 + 10*1.0
	assert like(0, 1) == np.log(1) - np.log(1) + 2*1.0 + 3.5*1.0
	assert like(1, 0) == np.log(2 - 2/3) - np.log(2/3) + 1.5*1.0 + 10.0*1.0
	assert like(1, 1) == np.log(2/3) - np.log(2 - 2/3) + 1.5*1.0 + 3.5*1.0

# Test simple (convex) optimizations
def test_optimize():
	ps = basic_ps()
	def like(p1, p2):
		return ps._optimize_cluster({'lig1':p1, 'lig2':p2}, 5)

	opt = {'lig1':1, 'lig2':0}
	assert like(0, 0)[1] == opt
	assert like(0, 1)[1] == opt
	assert like(1, 0)[1] == opt
	assert like(1, 1)[1] == opt

def test_max():
	ps = basic_ps()

	opt = {'lig1':1, 'lig2':0}
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt

def three_ligands():
	ligands = {'lig1': create_ligand('lig1', [0, 1], [-2, -1.5],
	                                 [{}, {(1, 23): 1.0}]),
			   'lig2': create_ligand('lig2', [0, 1], [-10, -3.5],
			                         [{(1, 23): 1.0}, {}]),
			   'lig3': create_ligand('lig3', [0, 1], [-1, -.5],
			                         [{}, {(1, 23): 0.5}])}
	stats = {'native': {'sb': DensityEstimate(domain = (0, 1)),
						'hbond': DensityEstimate(domain = (0, 1))},
			 'reference': {'sb': DensityEstimate(domain = (0, 1)),
			 			   'hbond': DensityEstimate(domain = (0, 1))}}
	stats['native']['sb'].x = np.linspace(0, 1, 100)
	stats['reference']['sb'].x = np.linspace(0, 1, 100)
	stats['native']['sb'].fx = np.linspace(0, 2, 100)
	stats['reference']['sb'].fx = np.linspace(2, 0, 100)
	features = {'sb': [1]}
	return PredictStructs(ligands, None, stats, features, 3, 1.0)

def test_max_three():
	ps = three_ligands()

	opt = {'lig1':1, 'lig2':0, 'lig3': 1}
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
	assert ps.max_posterior() == opt
