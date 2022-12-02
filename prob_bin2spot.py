from head_lib import *

def compute_threshold(tmap_adata):
	n_bin = tmap_adata.n_obs
	bin_value_list = []
	for r_index in range(n_bin):
		b=tmap_adata.X[r_index].copy()
		c_cell_index = b.argsort()[-1:]
		c_value = b[c_cell_index][0]
		bin_value_list.append(c_value)

	bin_value_list.sort()
	bin_sum = sum(bin_value_list)
	threshold = 0
	c_sum = 0
	for i in bin_value_list:
		c_sum+=i
		if c_sum/bin_sum >= 0.005:
			threshold = i
			break

	return threshold


def trans_porb_bin2spot(map_adata,spot_adata,featrue,stored_num):
	"""
	Transfrom the bin size cell_ID probability to spot's cell_ID proability
	Args:
		map_adata(AnnData): The mapping matrix in anndata format by tangram
		spot_adata(AnnData): The gene expression matrix for spot in anndata format.
		featrue (String): The featrue be assigned ("cell_ID")
		stored_num(Int) : The top probability mark ("cell_ID") to store for each spot.
	Retrun:
		index_adata,adata(AnnData) : The 2 anndatas to stored spot's most likely cell_ID and probability values
	"""
	tmap_adata = map_adata.T
	bin_df = tmap_adata.obs
	n_bin = tmap_adata.n_obs
	n_obs = spot_adata.n_obs
	n_var = stored_num

	row = np.zeros(n_var)
	col = np.arange(n_var)
	index_data = tmap_adata.X[0].argsort()[-n_var:] ####
	data = tmap_adata.X[0][index_data]
	X_matrix = sparse.coo_matrix((data,(row,col)),shape=(n_obs,n_var))
	X_matrix = X_matrix.tolil()
	index_matrix = sparse.coo_matrix((index_data,(row,col)),shape=(n_obs,n_var))
	index_matrix = index_matrix.tolil()
	new_obs = spot_adata.obs
	top_cell = np.arange(n_var)
	new_var = pd.DataFrame(data = top_cell,columns=['top_cell'])

	threshold = compute_threshold(tmap_adata)

	for r_index in range(n_bin):
		bin_ID = bin_df.loc[str(r_index),featrue]
		b=tmap_adata.X[r_index].copy()
		cell_index = b.argsort()[-n_var:]
		line = b[cell_index]
		index = new_obs.loc[new_obs[featrue]==bin_ID].index.astype(int)
		for i in index:
			index_matrix[i] = cell_index
			X_matrix[i] = line

	X_matrix = X_matrix.tocsr()
	index_matrix = index_matrix.tocsr()
	adata = ad.AnnData(X=X_matrix,obs=new_obs,var=new_var)
	index_adata = ad.AnnData(X=index_matrix,obs=new_obs,var=new_var)
	
	return index_adata,adata


def trans_porb_bin2spot_v2(map_adata,spot_adata,featrue,stored_num):
	tmap_adata = map_adata.T
	bin_df = tmap_adata.obs
	n_bin = tmap_adata.n_obs
	n_obs = spot_adata.n_obs
	n_var = stored_num

	row = np.zeros(2)
	col = np.arange(2)
	index_data = tmap_adata.X[0].argsort()[-n_var:][0] ####
	data = tmap_adata.X[0][index_data]
	in_va_data = np.array([index_data,data])

	X_matrix = sparse.coo_matrix((in_va_data,(row,col)),shape=(n_obs,2))
	X_matrix = X_matrix.tolil()
	new_obs = spot_adata.obs.copy()
	new_obs["cell_ID"] = -1
	top_cell = np.arange(2)
	new_var = pd.DataFrame(data = top_cell,columns=['top_cell_ID'])
	threshold = compute_threshold(tmap_adata)

	for r_index in range(n_bin):
		bin_ID = bin_df.loc[str(r_index),featrue]
		b=tmap_adata.X[r_index].copy()
		cell_index = b.argsort()[-n_var:][0]
		value = b[cell_index]
		index = new_obs.loc[new_obs[featrue]==bin_ID].index.astype(int)
		if value <= threshold:
			line = np.array([-1,0])
			bin_ID = -1
		else:
			line = np.array([cell_index,value])

		for i in index:
			X_matrix[i]=line
			new_obs.loc[i,["cell_bin","cell_ID"]] = [bin_ID,line[0]]

	X_matrix = X_matrix.tocsr()
	adata = ad.AnnData(X=X_matrix,obs=new_obs,var=new_var)
	
	return adata,threshold


def trans_porb_bin2spot_v3(map_adata,spot_adata,featrue):
	"""
	Transfrom the bin size cell_ID probability to spot's cell_ID proability
	Args:
		map_adata(AnnData): The mapping matrix in anndata format by tangram
		spot_adata(AnnData): The gene expression matrix for spot in anndata format.
		featrue (String): The featrue be assigned ("cell_ID")
	Retrun:
		index_adata,adata(AnnData) : The 2 anndatas to stored spot's most likely cell_ID and probability values
	"""
	tmap_adata = map_adata.T
	bin_df = tmap_adata.obs
	n_var = tmap_adata.n_vars
	n_bin = tmap_adata.n_obs
	n_obs = spot_adata.n_obs
	row = np.zeros(n_var)
	col = np.arange(n_var)
	data = tmap_adata.X[0]
	X_matrix = sparse.coo_matrix((data,(row,col)),shape=(n_obs,n_var))
	X_matrix = X_matrix.tolil()
	new_obs = spot_adata.obs
	top_cell = np.arange(n_var)
	new_var = pd.DataFrame(data = top_cell,columns=['top_cell'])

	for r_index in range(n_bin):
		bin_ID = bin_df.loc[str(r_index),featrue]
		b=tmap_adata.X[r_index].copy()
		index = new_obs.loc[new_obs[featrue]==bin_ID].index.astype(int)
		for i in index:
			X_matrix[i] = b

	X_matrix = X_matrix.tocsr()
	adata = ad.AnnData(X=X_matrix,obs=new_obs,var=new_var)
	
	return adata


def trans_porb_bin2spot_v4(map_adata,spot_adata,featrue,stored_num):
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
