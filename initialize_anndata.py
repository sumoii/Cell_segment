import sys
from head_lib import *
from us_log import *
from run_tangram import *
from prob_bin2spot import *

def namelist2numlist(name_list):
	"""
	Generate correspoding numerical list from a name list. The numerical id of a name is the order of the name appearing in the name list. The numerical id starts from 0. 
	Args：
		name_list(list): a list of string
	Return:
		num_list(list): 
	"""
	num_list = []
	unique_name_dict = {}
	unique_name_list = []
	order_id = 0
	for name in name_list:
		if name not in unique_name_dict.keys():
			unique_name_dict[name] = order_id
			unique_name_list.append(name)
			order_id = order_id + 1
		
		num_list.append(unique_name_dict[name])

	unique_name_num  = order_id  

	return num_list, unique_name_list, unique_name_num
	

def gem2df(gemfile):
	"""
	Read a gem file to dataframe and add a new column of 'x_y'.
	Args:
		gemfile(str): file name with path
	Return:
		df_gem(DataFrame): the datafrme
	"""
	df_gem = pd.read_csv(gemfile, delim_whitespace=True)
	df_gem = df_gem.sort_values(by=['x','y'])
	df_gem['x_y'] = df_gem['x'].map(str) + "_" + df_gem['y'].map(str)

	return df_gem


def is_in_area(x,range_str):
	"""
	Judge the x is in the range
	Args:
		x(int): A spot's coordinate
		range_str(str) : A given range
	Retturn:
		True or False
	"""
	start = int(range_str.split("_")[0])
	end = int(range_str.split("_")[1])
	data = int(x)
	if start <= data and data<= end :
		return True

	return False


def mark_cell_area(spot_matrix,bin_size):
	"""
	Marking cell_ID for the spot_matrix
	Args:
		spot_matrix(DataFrame): A dataframe from gem2df
		bin_size(Int) : A bin size
		x_blank,y_blank(Array) : A array of the blank(track) lines
	Retrun:
		area_cell_type(Dictionary) : A dictionary with cell area(bin)
	"""
	s_bin = bin_size - 1
	a_max = spot_matrix.max();a_min = spot_matrix.min()
	x_start = int(a_min[0]);x_end = int(a_max[0])
	y_start = int(a_min[1]);y_end = int(a_max[1])
	area_cell_type = {}
	#x_blank,y_blank = locate_blank_line(spot_matrix,x_start,x_end,y_start,y_end)
	cell_count = 0
	x_s = x_start;x_e = x_start;y_s = y_start;y_e = y_start
	spot_xlist = []
	spot_ylist = []
	spot_index = []
	spot_cell_ID =[]
	while x_e <= x_end:	
		x_e += s_bin
		if x_e >= x_end or x_end - x_e < (bin_size/4):
			x_e=x_end
		#x_e,x_region = delete_bin_trackline(x_s,x_e,x_blank)
		x_region = np.arange(x_s,x_e+1)
		while y_e <= y_end:
			y_e += s_bin
			if y_e >= y_end or y_end - y_e < (bin_size/4):
				y_e = y_end
			#y_e,y_region = delete_bin_trackline(y_s,y_e,y_blank)
			y_region = np.arange(y_s,y_e+1)
			for x in x_region:
				for y in y_region:
					spot_xlist.append(x)
					spot_ylist.append(y)
					row_index = str(x)+"_"+str(y)
					spot_index.append(row_index)
					spot_cell_ID.append(cell_count)
			cell_count+=1
			y_e += 1;y_s = y_e
		x_e += 1;x_s = x_e
		y_s = y_start;y_e = y_start

	spot_data = {"x":spot_xlist,"y":spot_ylist,"cell_bin":spot_cell_ID}

	new_df = pd.DataFrame(data=spot_data,index=spot_index)

	return new_df,x_start,x_end,y_start,y_end


def add_spot_to_matrix(new_df,old_df,spot_gene_matrix):
	new_df["is_add_spot"] = 1
	new_df.loc[old_df.index,"is_add_spot"] = 0
	new_df = new_df.reset_index()
	new_df = new_df.sort_values(by=["index"])
	new_df = new_df.reset_index(drop=True)
	new_spot_index = new_df.loc[new_df["is_add_spot"]==0].index
	spot_num = old_df.shape[0]
	new_spot_num = new_df.shape[0]
	row_index = np.zeros(spot_num)
	col_index = np.arange(spot_num)
	data = np.zeros(spot_num)
	new_spot_spot = sparse.coo_matrix((data,(row_index,col_index)),shape=(new_spot_num,spot_num))
	new_spot_spot = new_spot_spot.tolil()

	for i in new_spot_index:
		old_row = old_df.index.get_loc(new_df.loc[i,"index"])
		new_spot_spot[i,old_row] = 1
	
	new_spot_spot = new_spot_spot.tocsr()

	new_spot_gene_matrix = new_spot_spot * spot_gene_matrix
	
	return new_df,new_spot_gene_matrix


def integrate_mark(mark_df,spot_num,gene_num,X_matrix,feature):
	"""
	integrate the spot's feature mark to a new sparse matrix.
	Args:
		mark_df(DataFrame): A dataframe for spot's marks.
		spot_num,gene_num(Int) : The numbers of all spot and gene.
		X_matrix(Sparse Matrix): The sparse matrix for spot's gene expression
		feature(String) : The feature mark of spot
	Return:
		new_X_matrix(Sparse Matrix): The sparse matrix for feature mark's gene expression
		new_obs(DataFrame): A datafrme for feature
	"""
	cell_ID_list =list(set(mark_df[feature].tolist()))
	if -1 in cell_ID_list:
		cell_ID_list.remove(-1)
	cell_ID_num = len(cell_ID_list)
	cell_ID_df_data= {feature:cell_ID_list}
	new_obs = pd.DataFrame(data=cell_ID_df_data)
	cell_ID_data = np.ones(gene_num)
	cell_ID_row = np.zeros(gene_num)
	cell_ID_col = np.arange(gene_num)
	new_X_matrix = sparse.csr_matrix((cell_ID_data,(cell_ID_row,cell_ID_col)),shape=(cell_ID_num,gene_num))
	new_X_matrix = new_X_matrix.tolil()
	for i in range(cell_ID_num):
		cell_ID = cell_ID_list[i]
		index = mark_df.loc[mark_df[feature]==cell_ID].index
		index = np.array(index).astype(int)
		spot_index = np.zeros(spot_num)
		spot_index[index] = 1
		result = spot_index * X_matrix
		new_X_matrix[i]=result
	
	new_X_matrix = new_X_matrix.tocsr()
	
	return new_obs,new_X_matrix


def calculate_bin_density(bin_size,genes_value):
	return genes_value/(bin_size*bin_size)


def delete_null_bin(bin_adata,spot_df,feature,bin_size,threshold,model):
	"""
	Delete the bin with no or low density gene expression.
	Args:
		bin_adata(AnnData): A AnnData for bin transcript data.
		bin_size(int): the initialize bin_size.
		spot_df(DataFrame):
		feature(String): which feature be specified.
		threshold(Float): the limit threshold.
		model(Int): The model be specified.
	Return:
	"""
	bin_num = bin_adata.n_obs
	bin_df = bin_adata.obs
	var_df = bin_adata.var
	X_matrix = bin_adata.X
	new_bin_list = []
	bin_index = []
	null_bin_list = []
	if model == 0:
		for i in range(bin_num):
			bin_name = bin_df.iloc[i][feature]
			bin_genes = X_matrix[i].sum()
			bin_density = calculate_bin_density(bin_size,bin_genes)
			if bin_density >= threshold:
				bin_index.append(i)
				new_bin_list.append(bin_name)
			else:
				null_bin_list.append(bin_name)

	elif model == 1:
		for i in range(bin_num):
			bin_name = bin_df.iloc[i][feature]
			bin_genes = X_matrix[i].sum()
			if bin_genes > 0:
				bin_index.append(i)
				new_bin_list.append(bin_name)
			else:
				null_bin_list.append(bin_name)

	new_bin_num = len(new_bin_list)
	row_index = np.zeros(bin_num)
	col_index = np.arange(bin_num)
	data = np.zeros(bin_num)
	temp_X_matrix = sparse.coo_matrix((data,(row_index,col_index)),shape=(new_bin_num,bin_num))
	temp_X_matrix = temp_X_matrix.tolil()
	for i in range(new_bin_num):
		old_row = bin_index[i]
		temp_X_matrix[i,old_row] = 1

	temp_X_matrix = temp_X_matrix.tocsr()
	bin_data = {feature:new_bin_list}
	obs_bin_df = pd.DataFrame(data = bin_data)
	new_X_matrix = temp_X_matrix * X_matrix
	bin_anndata = ad.AnnData(X=new_X_matrix, obs= obs_bin_df ,var=var_df)

	null_bin_num = len(null_bin_list)
	for i in range(null_bin_num):
		null_bin_name = null_bin_list[i]
		index = spot_df.loc[spot_df[feature]==null_bin_name].index
		spot_df.loc[index,feature] = -1

	return bin_anndata,spot_df
	

def gem2anndata(df_gem,bin_size,threshold):
	"""
	Constructing anndata based on the information in Dataframe and marking cell lable.
	Args:
		df_gem(DataFrame): A dataframe from gem2df.
		bin_size(Int) : A int of bin size.
		threshold(Int) : A int of the bin gene density limit.
	Return:
		t_anndata(Anndata): Anndata where X is a matrix in CSR format, obs is the spot, and var is geneID.
	"""
	gene_list = df_gem['geneID'].tolist()
	spot_list = df_gem['x_y'].tolist()
	count_list = df_gem['MIDCounts'].tolist()
	gene_num_list, gene_unique_list, gene_num = namelist2numlist(gene_list)
	spot_num_list, spot_unique_list, spot_num = namelist2numlist(spot_list)
	row = np.array(spot_num_list)
	col = np.array(gene_num_list)
	data = np.array(count_list)
	X_matrix = sparse.coo_matrix((data, (row, col)), shape=(spot_num, gene_num))
	X_csr = X_matrix.tocsr()
	spot_xlist = [x.split("_")[0] for x in spot_unique_list]
	spot_ylist = [x.split("_")[1] for x in spot_unique_list]
	spot_data = {"x": spot_xlist, "y": spot_ylist}
	df_obs = pd.DataFrame(data=spot_data, index=spot_unique_list)
	df_obs['x'] = pd.to_numeric(df_obs['x'])
	df_obs['y'] = pd.to_numeric(df_obs['y'])
	new_df_obs,x_start,x_end,y_start,y_end= mark_cell_area(df_obs,bin_size)
	new_df_obs,new_spot_gene_matrix = add_spot_to_matrix(new_df_obs,df_obs,X_csr)
	new_spot_num = new_df_obs.shape[0]
	df_cell_obs,cell_csr = integrate_mark(new_df_obs, new_spot_num, gene_num, new_spot_gene_matrix, "cell_bin")
	df_var = pd.DataFrame(data=gene_unique_list, index=gene_unique_list, columns=['geneID'])
	bin_anndata = ad.AnnData(X=cell_csr, obs=df_cell_obs, var=df_var)
	bin_anndata,new_df_obs = delete_null_bin(bin_anndata, new_df_obs, "cell_bin", bin_size, threshold, 1)
	t_anndata = ad.AnnData(X=new_spot_gene_matrix, obs=new_df_obs, var=df_var)

	return x_start,x_end,y_start,y_end,t_anndata,bin_anndata


def initialize_mrf(x_start,x_end,y_start,y_end,spot_df):
	"""
	Initialize the vector of the spot
	Args:
		spot_df : The DataFrame of spot's coordinate data.

	Return:
		mrf_adata(AnnData) : A annData for spot's index ,which X sotred index.
		which can accelerate the subsequent processing steps --findnewbin.
	"""
	x_array = np.arange(x_start,x_end+1)
	y_array = np.arange(y_start,y_end+1)
	x_data = {"x":x_array}
	y_data = {"y":y_array}
	x_shape = len(x_array)
	y_shape = len(y_array)
	x_num = []
	y_num = []
	data_num = []
	df = spot_df.copy()
	spot_num = spot_df.shape[0]
	for i in range(spot_num):
		x = df.iloc[i]["x"]
		y = df.iloc[i]["y"]
		x = x - x_start
		y = y - y_start
		x_num.append(x)	
		y_num.append(y)
		data_num.append(i)

	spot_matrix = sparse.coo_matrix((data_num,(x_num,y_num)),shape=(x_shape,y_shape))
	spot_matrix = spot_matrix.tocsr()
	df_obs = pd.DataFrame(data=x_data)
	df_var = pd.DataFrame(data=y_data)
	mrf_adata = ad.AnnData(X=spot_matrix,obs=df_obs,var=df_var)

	return mrf_adata


def output_anndata(t_anndata, outfile):
	"""
	Output an anndata into a specified file.
	Args:
		t_anndata(Anndata): Anndata
		outfile(str): output filename with path information
	Return: None
	"""
	t_anndata.write(outfile)
