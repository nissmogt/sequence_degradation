#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import sys
import numpy as np
from analysis.analysis_pipeline import pipeline_replicates, process_dca
from msa.tools.check_length import check_length, check_nseq
import data.tools.pdb
from analysis.plots import (plot_average_ppv, multiple_plot_average_ppv, plot_avg_zscore, plot_neff_vs_zscore,
                            plot_avg_dist, plot_fraction_below_threshold, plot_ptp)

# sysid = sys.argv[1].strip(".fa")
# root = os.path.join("/scratch", "kmm5")
# dir_pdb = os.path.join("/scratch", pdb")
sysid = "1cc8A"
root = os.path.join("tests", "single")
dir_pdb = "pdb"
print(f"System ID: {sysid}\nPath: {root}")
dir_sys = os.path.join(root, "systems", sysid)
dir_out = os.path.join(dir_sys, "results")
msa_file = os.path.join(dir_sys, f"{sysid}_filtered_25.fasta")

neff_file = os.path.join(dir_sys, "neff_array.npy")
neff = np.load(neff_file)

seq_l = check_length(msa_file)
nseq = check_nseq(msa_file)
threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]

df_pdb = data.tools.pdb.pipeline_pdb(sysid, dir_pdb)
df_z = process_dca(root=root, _sysid=sysid, _df_pdb=df_pdb, _nseqs=nseq, _neff=neff, _rep=1, zcalc=True,
                   shift=4, replicates=False)
