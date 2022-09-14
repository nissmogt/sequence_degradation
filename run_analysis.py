#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import sys
import numpy as np
from analysis.analysis_pipeline import pipeline_replicates
from msa.tools.check_length import check_length
from analysis.plots import (plot_average_ppv, multiple_plot_average_ppv, plot_avg_zscore, plot_neff_vs_zscore,
                            plot_avg_dist, plot_fraction_below_threshold, plot_ptp)

sysid = sys.argv[1].strip(".fa")
#sysid = "1a3aA"
dca_dir = os.path.join("/scratch", "kmm5")
print(f"System ID: {sysid}\nPath: {dca_dir}")
msa_file = os.path.join(dca_dir, "systems", sysid, f"{sysid}_filtered_25.fasta")
dir_out = os.path.join(os.path.dirname(msa_file), "results")
neff_file = os.path.join(os.path.dirname(msa_file), f"replicates/neff_array.npy")
msa = os.path.join(dca_dir, "systems", sysid, f"{sysid}_filtered_25.fasta")
seq_l = check_length(msa_file)
threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]

ppv, pred_rank, top_z, top_dist = pipeline_replicates(dca_dir, sysid, seq_l, threshold_values, npairs=0, zfilter=True,
                                                      plots=True, passthrough=False)
neff = np.load(neff_file)
neff = neff[:, :-1]
avg_neff = np.mean(neff[:60, :], axis=0)
avg_dist = np.mean(top_dist[:60, :, :], axis=0)
plot_fraction_below_threshold(avg_dist, avg_neff, seq_l, sysid, dir_out, extra_text="", save=True)
plot_avg_zscore(top_z[:60, :, :], avg_neff, sysid, seq_l, dir_out, extra_text="", save=True)
plot_neff_vs_zscore(top_z[:60, :, :], avg_neff, sysid, seq_l, dir_out, extra_text="", save=True)
plot_ptp(top_z[:60, :, :], avg_neff, sysid, seq_l, dir_out, extra_text="", save=True)
plot_avg_dist(top_dist[:60, :, :], 8, avg_neff, sysid, seq_l, dir_out, extra_text="", save=True)
multiple_plot_average_ppv(ppv[:, :60, :], neff[:60, :], sysid, seq_l, threshold_values, dir_out,
                          norm=1, extra_text="8A", save=True)
