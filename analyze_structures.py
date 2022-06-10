import glob
import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd
from Bio.PDB import *
from distance_pdb import distance_matrix, import_pdb_distance_matrix

# from prody import *
# 1) Load pdb

# folder = "structures\\p300_CBP_TAZ1\\"
# for i, s in enumerate(glob.glob(folder+'*.ent')):
#     ch = s.split("\\")[-1].strip(".ent")[-1]
#     df, chn = distance_matrix(s, all_atom=True, chain=ch)

# average distance between same indices
folder = "contact_maps\\"
cutoff = 8
df = pd.read_csv(folder + "concatenated_taz1_distance_matrix.txt", delimiter="\t")
df_filter = df[df["d"] <= cutoff]
sub_df = df_filter.iloc[:, :3]
average_d = sub_df.groupby(["i", "j"]).mean().reset_index()
header = ["i", "j", "d"]
average_d.to_csv("taz1_avg_distance_matrix_8A.txt", sep='\t', index=False, header=header, float_format='%.5f')


# outfile = folder+"concatenated_taz1_distance_matrix.txt"
# f = open(folder+"atom_distance_matrix_1l3e_B.ent.txt", "r")
# header = f.readline()
# f.close()
# l = []
# l.append(header)
# count = 0
# for i, s in enumerate(glob.glob(folder+'*.txt')):
#     with open(s, 'r') as f:
#         f.readline()
#         for j, line in enumerate(f.readlines()):
#             l.append(line)
#         count += j
#         print(count)
# with open(outfile, "w", newline="") as f:
#     for i, line in enumerate(l):
#         f.writelines(line)
#     print(i-1)
# print(i-1-count-8)

# df_pdb = import_pdb_distance_matrix(s.split("\\")[-1], all_atom=True)
# df_filtered = df_pdb[df_pdb["dist"] <= cutoff]


# n_nmr = len(structure)
# u = parsePDB(folder+fpdb, subset='ca', chain='A')
# ensemble = Ensemble('Taz1 NMR ensemble')
# ensemble.setCoords(u.getCoords())
# ensemble.addCoordset(u.getCoordsets())
# ensemble.iterpose()
# u_copy = u.copy()
# u_copy.delCoordset(range(u_copy.numCoordsets()))
# u_copy.addCoordset(ensemble.getCoordsets())
#
# pca = PCA('Taz1')
# pca.buildCovariance(ensemble)
# pca.calcModes()
# repr(pca)
