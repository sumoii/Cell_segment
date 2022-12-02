from head_lib import *

def trans_porb_bin2spot(map_adata,spot_adata,featrue,stored_num):
	"""
	Transfrom the bin size cell_ID probability to spot's cell_ID proability
	Args:
		map_adata(AnnData): The mapping matrix in anndata format by tangram.
		spot_adata(AnnData): The gene expression matrix for spot in anndata format.
		featrue (String): The featrue be assigned.
	Retrun:
		adata(AnnData): The anndatas stored spot's most  front(stored_num )cell_ID's probability values.
	"""
	tmap_adata = map_adata.T
	bin_df = tmap_adata.obs
	n_bin = tmap_adata.n_obs
	n_obs = spot_adata.n_obs
	n_var = stored_num
	row = np.zeros(n_var)
	col = np.arange(n_var)
	data = np.zeros(n_var)
	X_matrix = sparse.coo_matrix((data,(row,col)),shape=(n_obs,n_var))
	X_matrix = X_matrix.tolil()
	new_obs = spot_adata.obs
	top_cell = np.arange(n_var)
	new_var = pd.DataFrame(data = top_cell,columns=['top_cell'])

	for r_index in range(n_bin):
		bin_ID = bin_df.loc[str(r_index),featrue]
		b=tmap_adata.X[r_index].copy()
		line = b[0:n_var]
		index = new_obs.loc[new_obs[featrue]==bin_ID].index.astype(int)
		for i in index:
			X_matrix[i] = line

	X_matrix = X_matrix.tocsr()
	adata = ad.AnnData(X=X_matrix,obs=new_obs,var=new_var)
	
	return adata
