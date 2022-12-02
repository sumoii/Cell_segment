from head_lib import *
from us_log import *
import scanpy as sc
import tangram as tg

def run_tangram(st_adata,sc_adata,group_name):
	"""
	Generate the mapping matrix by tangram for single cell data and spatial transcript data
	Args:
		st_adata(AnnData) : The anndata of spatial transcript data 
		sc_adata(AnnData) : The  of spatial transcript data
		group_name(String) : The group name in single cell data for chooese traning genes in preprocessing
	Return:
		ad_map(AnnData) : The mapping matrix in anndata format
	"""
	sc.pp.normalize_total(st_adata)
	sc.pp.log1p(st_adata)
	sc.tl.rank_genes_groups(sc_adata, groupby=group_name, use_raw=False)
	markers_df = pd.DataFrame(sc_adata.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
	markers  =  list(np.unique(markers_df.melt().value.values))
	tg.pp_adatas(sc_adata,st_adata, genes=markers)

	ad_map = tg.map_cells_to_space(sc_adata, st_adata, mode="cells", density_prior='rna_count_based', num_epochs=500, device='cpu')

	return ad_map
