import queue
import sys
import pandas as pd
import anndata as ad

def is_in_area(x,y,x_end,y_end):
	if 0<=x and x<=x_end and 0<=y and y<=y_end:
		return True
	
	return False


def BFS(cell_ID_df,X_matrix,x_start,x_end,y_start,y_end):
	cell_bin_num = 0
	boundary_bin_id = []
	new_x_end = x_end - x_start
	new_y_end = y_end - y_start
	x_list = [-1,0,1,-1,1,-1,0,1]
	y_list = [1,1,1,0,0,-1,-1,-1]
	nextstartspot = queue.Queue()
	nextstartspot.put("0_0")
	q = queue.Queue()
	bin_num_dict = {}
	while not nextstartspot.empty():
		q_start = nextstartspot.get()
		q.put(q_start)
		while not q.empty():
			u = q.get()
			x = int(u.split("_")[0])
			y = int(u.split("_")[1])
			now_df_index = int(X_matrix[x,y])
			if cell_ID_df.loc[now_df_index,"visited"] == 1:
				continue
			cell_ID_df.loc[now_df_index,"visited"] = 1
			cell_ID_df.loc[now_df_index,"new_cell_bin"] = cell_bin_num
			if cell_ID_df.loc[now_df_index,"new_cell_ID"] == 0:
				boundary_bin_id.append(cell_bin_num)
			for j in range(8):
				new_x = x+x_list[j]
				new_y = y+y_list[j]
				if is_in_area(new_x,new_y,new_x_end,new_y_end):
					v = str(new_x)+"_"+str(new_y)
					new_df_index = int(X_matrix[new_x,new_y])
					if cell_ID_df.loc[new_df_index,"visited"] == 1:
						continue
					if cell_ID_df.iloc[now_df_index]["new_cell_ID"]==cell_ID_df.iloc[new_df_index]["new_cell_ID"]:
						q.put(v)
					else:
						nextstartspot.put(v)
		cell_bin_num+=1

	return cell_ID_df,boundary_bin_id


def find_new_cell_bin(cell_ID_df,mrf_2d,x_start,x_end,y_start,y_end):
	X_matrix = mrf_2d.X
	cell_ID_df["visited"] = 0
	cell_ID_df["new_cell_bin"] = -1
	cell_ID_df,boundary_bin_id = BFS(cell_ID_df,X_matrix,x_start,x_end,y_start,y_end)

	for i in boundary_bin_id:
		index = cell_ID_df[cell_ID_df["new_cell_bin"]==i].index
		cell_ID_df.loc[index,"new_cell_bin"] = -1

	bin_list = list(set(cell_ID_df["new_cell_bin"].tolist()))
	bin_list.remove(-1)

	for j in range(len(bin_list)):
		bin_num = bin_list[j]
		index = cell_ID_df.loc[cell_ID_df["new_cell_bin"]==bin_num].index
		cell_ID_df.loc[index,"new_cell_bin"] = j

	return cell_ID_df


if __name__=="__main__":
	mrf_2d = ad.read_h5ad(sys.argv[1])
	region = sys.argv[2]
	cell_ID_df = pd.read_csv(sys.argv[3],delim_whitespace=True)
	prefix = sys.argv[4]
	X_matrix = mrf_2d.X
	x_start = int(region.split("_")[0])
	x_end = int(region.split("_")[1])
	y_start = int(region.split("_")[2])
	y_end = int(region.split("_")[3])
	cell_ID_df["visited"] = 0
	cell_ID_df["new_cell_bin"] = -1
	cell_ID_df = BFS(cell_ID_df,X_matrix,x_start,x_end,y_start,y_end)

	bin_list = list(set(cell_ID_df["new_cell_bin"].tolist()))
	for j in range(len(bin_list)):
		bin_num = bin_list[j]
		index = cell_ID_df.loc[cell_ID_df["new_cell_bin"]==bin_num].index
		cell_ID_df.loc[index,"new_cell_bin"] = j

	cell_ID_df.to_csv(prefix,sep="\t")
