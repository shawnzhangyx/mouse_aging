#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='Run NMF using sklearn.')
parser.add_argument('--normH', type=str, dest="normH", help='coefficient matrix H')
parser.add_argument('--statH', type=str, dest="statH", help='input statH matrix')
parser.add_argument('-p', '--perplexity', type=int, dest="perplexity", default=0.3, help='an int for perplexity')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')

args = parser.parse_args()

from os.path import dirname, abspath, join
from warnings import warn

import numpy as np
import scipy as sp
from scipy import io

from sklearn.metrics import silhouette_samples, silhouette_score, pairwise_distances
import umap

from time import perf_counter as pc

try:
    from matplotlib.pyplot import savefig, imshow, set_cmap
except ImportError as exc:
    warn("Matplotlib must be installed.")

def run():
	""" Run standard NMF on rank """
	start_time = pc()
	""" init input files """
	normHf = args.normH
	statHf = args.statH
	outPrefix = args.outPrefix
	perplexity = args.perplexity
	inF = read_files(normHf, statHf)
	normH = inF['normH']
	o_stat_H = inF['statH']
	rank = inF['rank']
	print("draw silhouette & tsne plot")
	X_dist = cal_pairwise_pearson(normH)
	X_transformed = cal_umap(X_dist)
	np.savetxt('.'.join([outPrefix, "umap.xy"]), X_transformed, fmt= "%g", delimiter="\t")
	
	end_time = pc()
	print('Used (secs): ', end_time - start_time)

def read_files(normH, statH):
	normH = np.loadtxt(normH)
	statH = np.genfromtxt(statH, dtype=None, names=True)
	rank = normH.shape[0]
	return {'normH':normH, 'statH':statH, 'rank':rank}
	
def cal_pairwise_dist(normH):
	X = normH.T
	X_pairwise_dist = pairwise_distances(X)
	return X_pairwise_dist

def cal_pairwise_pearson(normH):
	X = normH.T
	X_corrcoef = np.corrcoef(X)
	X_corrcoef_dist = np.sqrt(2*(1-X_corrcoef))
	return X_corrcoef_dist

def cal_umap(X):
	embedding = umap.UMAP(metric='correlation')
	#embedding = umap.UMAP(n_neighbors=p, min_dist=0.1, metric='correlation')
	X_transformed = embedding.fit_transform(X)
	return X_transformed

if __name__ == "__main__":
	"""Draw U-MAP plot for coefficient matrix"""
	run()
