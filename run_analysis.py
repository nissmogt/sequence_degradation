#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import sys
import numpy as np
from analysis.analysis_pipeline import pipeline_replicates, process_dca
from msa.tools.check_length import check_length, check_nseq
import data.tools.pdb
from analysis.plots import (plot_average_ppv, multiple_plot_average_ppv, plot_avg_zscore, plot_neff_vs_zscore,
                            plot_avg_dist, plot_fraction_below_threshold, plot_ptp)

sysid = sys.argv[1].strip(".fa")
root = os.path.join("/scratch", "kmm5", "single")
dir_pdb = os.path.join("/scratch", "kmm5", "pdb")
# sysid = "1cc8A"
# root = os.path.join("tests", "single")
# dir_pdb = "pdb"
print(f"System ID: {sysid}\nPath: {root}")
dir_sys = os.path.join(root, "systems", sysid)
dir_out = os.path.join(dir_sys, "results")
msa_file = os.path.join(dir_sys, f"{sysid}_filtered_25.fasta")

neff_file = os.path.join(dir_sys, "neff_array.npy")
neff = np.load(neff_file)

seq_l = check_length(msa_file)
nseq = check_nseq(msa_file)
thresholds_list = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]

df_pdb = data.tools.pdb.pipeline_pdb(sysid, dir_pdb)
if sysid == "1cc8A":
    df_dca = process_dca(root=root, _sysid=sysid, _df_pdb=df_pdb, _nseqs=nseq, _neff=neff, _rep=1, zcalc=True,
                         shift=4, replicates=False)
else:
    df_dca = process_dca(root=root, _sysid=sysid, _df_pdb=df_pdb, _nseqs=nseq, _neff=neff, _rep=1, zcalc=True,
                         shift=0, replicates=False)

n_degraded_seqs = nseq
_ncols = seq_l
dir_contact_map = os.path.join(dir_sys, "images")
if not os.path.exists(dir_contact_map):
    os.makedirs(dir_contact_map)
prefix = ""
plots = True
zfilter = True
dca_score = "DIAPC"
distance_cutoff = 8
map_idx = 0
import analysis.plots
from analysis.validation import calculate_ppv

if plots:
    # plot_score_distribution(dca_score, n_degraded_seqs, df_dca, neff, ncols, dir_contact_map)
    if zfilter:
        print("zfilter")
        analysis.plots.plot_dist_distribution(df_dca[:750], n_degraded_seqs,
                                              _ncols, dir_contact_map, extra_text=prefix)
        analysis.plots.plot_score_distribution("zscore", n_degraded_seqs, df_dca, neff, _ncols,
                                               dir_contact_map, extra_text=prefix)

pos_pred_list = np.zeros(len(thresholds_list))
pair_rank_array = np.zeros_like(pos_pred_list)
if plots:
    # plot contact map without filter
    analysis.plots.compare_contact_map(df_dca, df_dca[:10], df_pdb, sysid, distance_cutoff,
                                       n_degraded_seqs, 0, _ncols, map_idx,
                                       dir_contact_map, extra_text="top_10")
for k in range(len(thresholds_list)):
    threshold_value = thresholds_list[k]
    if zfilter:
        df_filtered = df_dca[df_dca["zscore"] >= threshold_value]
        topn_z_array = df_dca["zscore"][:10].to_numpy()
        topn_dist_array = df_dca["d"][:10].to_numpy()
        if plots:
            print(f"{threshold_value}")
            analysis.plots.plot_top_zscore(df_dca, 10, n_degraded_seqs, _ncols,
                                           dir_contact_map, extra_text=prefix)
            # calculate average z-score and std for the top 10 di pairs
            analysis.plots.compare_contact_map(df_dca, df_filtered, df_pdb, sysid, distance_cutoff,
                                               n_degraded_seqs, threshold_value, _ncols, map_idx,
                                               dir_contact_map, extra_text=prefix)
    else:
        df_filtered = df_dca[:threshold_value]

    if len(df_filtered) > 0:
        ppv, pair_rank = calculate_ppv(df_filtered, distance_cutoff)
        if plots:
            analysis.plots.plot_ppv(zfilter, ppv, neff, sysid, _ncols,
                                    threshold_value, dir_contact_map, extra_text=prefix)
            print("plotting PPV")
        pos_pred_list[k] = ppv[-1]
        pair_rank_array[k] = pair_rank
    else:
        pos_pred_list[k] = 0

output_ppv = os.path.join(dir_out, f"{sysid}_ppv_zscore_{dca_score}_100reps.npy")
output_pair_rank = os.path.join(dir_out, f"{sysid}_ppv_pairrank_z_{dca_score}_reps100.npy")
output_z = os.path.join(dir_out, f"avg_z_{sysid}_100reps.npy")
output_dist = os.path.join(dir_out, f"avg_dist_{sysid}_100reps.npy")
np.save(output_ppv, pos_pred_list)
np.save(output_pair_rank, pair_rank_array)
np.save(output_z, topn_z_array)
np.save(output_dist, topn_dist_array)
