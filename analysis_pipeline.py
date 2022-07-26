import os
import numpy as np
import data.tools.pdb as dpdb
import analysis.dca_object as dca
import analysis.plots
from analysis.validation import calculate_ppv

# sysid = sys.argv[1].strip(".fa")
sysid = "1a3aA"
root = os.getcwd()
nseqs = 5165
raw_dca = f"DI_{sysid}_n{nseqs}.txt"
dir_dca = os.path.join(root, "systems", sysid, "replicates", "sub0", "mf", "pc0.2")
dca_in = os.path.join(dir_dca, raw_dca)

d = dca.DirectCoupling()
d._sysid = sysid
d._score = "diapc"
neff_array = np.load(os.path.join(root, "systems", sysid, "replicates", "neff_array.npy"))
n_effective = neff_array[0][0]
out_dca = os.path.join(dir_dca, f"{sysid}_neff{n_effective}_all.txt")

df = d.load_raw(dca_in)
df_rank = d.rank_hamming(df, distance=5)
df_pdb = dpdb.pipeline_pdb(sysid, "pdb")
df_dca = df_rank.merge(df_pdb, how='inner', on=['i', 'j'])

zscore_reference = np.load(os.path.join(root, "assets", f"monomer_DIAPC_scores.npy"))
df_z = d.zscore(df_rank, zscore_reference)

# d.savetocsv(df_z, out_dca)
analysis.plots.contact_map_single(df_z, monomer=df_pdb, n=10, x="test", distance_cutoff=8)
# TODO: Plot ppv vs zscore for every neff/l value
ppv = calculate_ppv(df_z, 6)
