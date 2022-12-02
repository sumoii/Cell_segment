import sys
import pydensecrf.densecrf as dcrf
import numpy as np
import anndata as ad
import pandas as pd


def create_pairwise_bilateral_2d(sx,sy,sr,img,W,H):
	"""
	Constructing the feature space of gene expression.
	Args:
		sx,sy,sr(float): sx,sy,sr the scaling value of feature space.
		img(array): The spot's genes counts.
		W,H(int): The width and height of the slice.
	Return:
		feats(array): The feature space array.
	"""
	feat_size = 3
	feats = np.zeros((feat_size,W,H),dtype=np.float32)
	for i in range(W):
		for j in range(H):
			feats[0, i, j] = i / sx
			feats[1, i, j] = j / sy
			feats[2, i, j] = img[i*H+j] / sr

	return feats.reshape([feat_size,-1])


def read_spot_Mcounts(spot_gene_ad):
	"""
	Read the spot gene count to a array.
	Args:
		spot_gene_ad(Anndata): The spot's gene expression matrix
	Return:
		img(array):The spot's gene count array.
	"""
	X_matrix = spot_gene_ad.X
	n_spot = spot_gene_ad.n_obs
	img = np.full(n_spot,0,dtype='int32')
	for i in range(n_spot):
		spot_sum = int(X_matrix[i].sum())
		img[i]=spot_sum

	return img


def read_spot_UnaryEnergy(value_ad,n_lables,spot_num):
	"""
	Read the spot's cell_ID probability matrix to two-dimensional array.
	Args:
		value_ad(Anndata): the spot's cell_ID probability matrix.
		n_lables(int) : the lables num is the cell num add one background.
		spot_num(int) : the slice's spot num.
	Return:
		U(arrat): two-dimensional array stored spot's cell_ID prob.
	"""
	U = np.full((n_lables,spot_num),0,dtype='float32')
	spot_num = value_ad.n_obs
	X_matrix = value_ad.X
	for i in range(spot_num):
		if X_matrix[i].sum()==0:
			back = -np.log(1)
		else:
			back = -np.log(0)
		U[0,i] = back
		for j in range(1,n_lables):
			U[j,i] = -np.log(X_matrix[i,j-1])

	return U


def run_densecrf(U, img, n_lables, W, H, w1, w2):
	"""
	Running Densecrf.
	Args:
		U(array): as above.
		img(array): as above.
		n_lables(int): as above.
		W,H(int): as above.
		w1,w2(int): as above.
	Return:
		MAP(array):two-dimensional array for spot's cell_ID prob. 
	"""
	d = dcrf.DenseCRF2D(H, W, n_lables)
	d.setUnaryEnergy(U)
	d.addPairwiseGaussian(sxy = (3,3), compat=w1, kernel=dcrf.DIAG_KERNEL, normalization=dcrf.NORMALIZE_SYMMETRIC)
	feats = create_pairwise_bilateral_2d(60,60,20, img, W, H)
	d.addPairwiseEnergy(feats, compat=w2, kernel=dcrf.DIAG_KERNEL, normalization=dcrf.NORMALIZE_SYMMETRIC)
	Q = d.inference(10)
	MAP = np.argmax(Q, axis=0)

	return MAP

