import pytest
import os
import sys
import numpy as np

from score.density_estimate import DensityEstimate

def test_average():
	de1 = DensityEstimate(points = 3, domain = (0, 2))
	de1.n_samples = 5
	de1.x = np.array([0, 1, 2])
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2))
	de2.n_samples = 15
	de2.x = np.array([0, 1, 2])
	de2.fx = np.array([1, 1, 1])

	avg = de1._average(de2)
	assert np.all(avg.fx == [0.75, 0.75, 0.75])

def test_average_zero():
	de1 = DensityEstimate(points = 3, domain = (0, 2))
	de1.n_samples = 0
	de1.x = np.array([0, 1, 2])
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2))
	de2.n_samples = 15
	de2.x = np.array([0, 1, 2])
	de2.fx = np.array([1, 1, 1])

	avg = de1._average(de2)

	assert np.all(avg.fx == [1, 1, 1])

def test_merge():
	de1 = DensityEstimate(points = 3, domain = (0, 2))
	de1.n_samples = 5
	de1.x = np.array([0, 1, 2])
	de1.fx = np.array([0, 0, 0])

	de2 = DensityEstimate(points = 3, domain = (0, 2))
	de2.n_samples = 15
	de2.x = np.array([0, 1, 2])
	de2.fx = np.array([1, 1, 1])

	de3 = DensityEstimate(points = 3, domain = (0, 2))
	de3.n_samples = 0
	de3.x = np.array([0, 1, 2])
	de3.fx = np.array([2, 2, 2])

	merged = DensityEstimate.merge([de1, de2, de3])

	assert np.all(merged.fx == [0.75, 0.75, 0.75])
