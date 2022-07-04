import os
import numpy as np
import sys
import analysis.dca_object as dca

# sysid = sys.argv[1].strip(".fa")
sysid = "1a3aA"
root = os.getcwd()
nseqs = 5165
raw_dca = f"DI_{sysid}_n{nseqs}.txt"
dir_dca = os.path.join(root, "systems", sysid, "replicates", "sub0", "mf", "pc0.2")

d = dca.DirectCoupling()
d._sysid = sysid
d._score = "diapc"
neff_array = np.load(os.path.join(root, "systems", sysid, "replicates", "neff_array.npy"))
n_effective = neff_array[0][0]
dca_in = os.path.join(dir_dca, raw_dca)

d.load_raw(dca_in)
d.rank_hamming()
df_dca = d._ranked_dataframe
dir_dca_in = os.path.dirname(dca_in)
outfile = os.path.join(dir_dca_in, f"{sysid}_neff{n_effective}_pc0.2_all.txt")

zscore_reference = np.load(os.path.join(root, f"monomer_{'diapc'.upper()}_scores.npy"))
# TODO: add reference file
df_z = d.zscore(df_dca, zscore_reference)
