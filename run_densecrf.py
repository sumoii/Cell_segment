import sys
import pydensecrf.densecrf as dcrf
import numpy as np
import anndata as ad
import pandas as pd


def create_pairwise_bilateral_2d(sx,sy,sr,img,W,H):
	feat_size = 3
	feats = np.zeros((feat_size,W,H),dtype=np.float32)
	for i in range(W):
		for j in range(H):
			feats[0, i, j] = i / sx
			feats[1, i, j] = j / sy
			feats[2, i, j] = img[i*H+j] / sr

	return feats.reshape([feat_size,-1])


def read_spot_Mcounts(spot_gene_ad):
	X_matrix = spot_gene_ad.X
	n_spot = spot_gene_ad.n_obs
	img = np.full(n_spot,0,dtype='int32')
	for i in range(n_spot):
		spot_sum = int(X_matrix[i].sum())
		img[i]=spot_sum

	return img


def read_spot_UnaryEnergy(X_df,n_lables,spot_num):
	U = np.full((n_lables,spot_num),0,dtype='float32')
	f = open(X_df,"r")
	index = 0
	line = f.readline()
	while line:
		line_list = line.split("\t")
		for i in range(n_lables):
			U[i,index] = line_list[i]
		line = f.readline()
		index+=1
	f.close()

	return U


def read_spot_UnaryEnergy_v2(value_ad,n_lables,spot_num):
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
	d = dcrf.DenseCRF2D(H, W, n_lables)
	d.setUnaryEnergy(U)
	d.addPairwiseGaussian(sxy = (3,3), compat=w1, kernel=dcrf.DIAG_KERNEL, normalization=dcrf.NORMALIZE_SYMMETRIC)
	feats = create_pairwise_bilateral_2d(60,60,20, img, W, H)
	d.addPairwiseEnergy(feats, compat=w2, kernel=dcrf.DIAG_KERNEL, normalization=dcrf.NORMALIZE_SYMMETRIC)
	Q = d.inference(10)
	MAP = np.argmax(Q, axis=0)

	return MAP

