import sys,argparse
from initialize_anndata import *
from run_densecrf import *
from findnewcellbin import *

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Cell segment (V1.0)')
	parser.add_argument('-g', required=True, dest='in_gem', type=str, help='the gem of stero-seq')
	parser.add_argument('-s', required=True, dest='sc_file', type=str, help='the h5ad of single cell')
	parser.add_argument('-b', required=True, dest='bin_size', type=int, help='the initialize size of a bin (spot number = n*n)')
	parser.add_argument('-d', required=True, dest='density', type=float, help='the gene density limit of a bin')
	parser.add_argument('-w1', required=True, dest='weight1', type=float, help='the weight of first Gaussian kernel')
	parser.add_argument('-w2', required=True, dest='weight2', type=float, help='the weight of second Gaussian kernel')
	parser.add_argument('-o', required=True, dest='out_file', type=str, help='the initialize output file prefix')
	parser.add_argument('-p', required=True, dest='prefix', type=str, help='the prefix of log file')
	parser.add_argument('-n', required=True, dest='t_time', type=int, help='the times of running tangram')
	parser.add_argument('-f', required=True, dest='feature', type=str, help='the feature input of tangram')

	args = parser.parse_args()	
	in_gem = args.in_gem
	sc_file = args.sc_file
	bin_size = args.bin_size
	threshold = args.density
	w1 = args.weight1
	w2 = args.weight2
	out_file = args.out_file
	prefix = args.prefix
	ct_time = args.t_time
	feature = args.feature
	sc_adata = ad.read_h5ad(sc_file)

	log_setting(prefix)
	logger.info("Start to initialize stereo-seq data to anndata")
	df_gem = gem2df(in_gem)
	x_start, x_end, y_start, y_end, spot_adata, bin_adata = gem2anndata(df_gem,bin_size,threshold)
	spot_df = spot_adata.obs
	spot_df.index = spot_df.index.astype(int)
	spot_adata.write(out_file+"_spot.h5ad")
	bin_adata.write(out_file+"_bin.h5ad")
	logger.info("===Finish initialize anndata")

	logger.info("Start to initialize cell_ID probability for each bin")
	mrf_adata = initialize_mrf(x_start,x_end,y_start,y_end,spot_df)
	map_adata = run_tangram(bin_adata,sc_adata,feature)
	mrf_adata.write(out_file+"_mrf.h5ad")
	map_adata.write(out_file+"_map.ha5d")
	logger.info("===Finish initialize cell_ID")

	logger.info("Start to cell segement")
	time = 1
	W = x_end-x_start+1
	H = y_end-y_start+1
	n_lables = sc_adata.n_obs+1
	gene_df = spot_adata.var
	spot_num = spot_adata.n_obs
	gene_num = spot_adata.n_vars
	X_matrix = spot_adata.X
	img = read_spot_Mcounts(spot_adata)
	value_adata = trans_porb_bin2spot(map_adata,spot_adata,"cell_bin",n_lables-1)
	U = read_spot_UnaryEnergy(value_adata,n_lables,spot_num)
	MAP = run_densecrf(U,img,n_lables,W,H,w1,w2)
	spot_df["new_cell_ID"] = MAP
	spot_df = find_new_cell_bin(spot_df,mrf_adata,x_start,x_end,y_start,y_end)
	spot_df.to_csv(str(time)+"_change.gem",sep="\t")
	time+=1
	while time<= ct_time:
		new_obs,new_X_matrix = integrate_mark(spot_df,spot_num,gene_num,X_matrix,"new_cell_bin")
		new_bin_adata = ad.AnnData(X=new_X_matrix,obs=new_obs,var=gene_df)
		new_bin_adata,spot_df = delete_null_bin(new_bin_adata,spot_df,"new_cell_bin",bin_size,threshold,1)
		spot_adata.obs = spot_df
		map_adata = run_tangram(new_bin_adata,sc_adata,feature)
		value_adata = trans_porb_bin2spot(map_adata,spot_adata,"new_cell_bin",n_lables-1)
		U = read_spot_UnaryEnergy(value_adata,n_lables,spot_num)
		MAP = run_densecrf(U,img,n_lables,W,H,w1,w2)
		spot_df["new_cell_ID"] = MAP
		spot_df = find_new_cell_bin(spot_df,mrf_adata,x_start,x_end,y_start,y_end)
		spot_df.to_csv(str(time)+"_change.gem",sep="\t")
		time+=1
	logger.info("===Finish cell segement")
