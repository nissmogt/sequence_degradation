import pandas as pd
import numpy as np
import networkx
# Construct an undirected graph where the edge is delta_td and the vertex is the residue index.

df = pd.read_csv("7lvs_1L8C_difference.csv", header=0)
delta_td = df.delta_td.to_numpy()
neighbor_list = np.array([df.i.to_numpy(), df.j.to_numpy()]).transpose()
x = [neighbor_list, delta_td]

# Centrality: indicators of centrality assign numbers or rankings to vertices within a graph
# corresponding to their network position

