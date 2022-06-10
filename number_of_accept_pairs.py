import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

plt.style.use('paper.mplstyle')


def check_length(filein):
    from Bio import AlignIO
    align = AlignIO.read(open(filein), "fasta")
    return align.get_alignment_length()


def process(filein):
    ppv_array = np.load(filein)
    avg_ppv = np.mean(ppv_array, axis=1)
    std_ppv = np.std(ppv_array, axis=1)
    n_thresholds, n_models = avg_ppv.shape
    x = np.zeros((n_thresholds, n_models))
    return ppv_array, avg_ppv, std_ppv, x


def plot_nap_tp(data_1, data_2, err, str_neff_l, thresh_vals, n_models, image_dir, dcutoff):
    nmodels, n_thresholds = data_1.shape
    cmap = plt.get_cmap("nipy_spectral")
    colors = [cmap(i) for i in np.linspace(0, 1, n_models + 1)]
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(5, 7))
    ax[0].set_ylabel('$\overline{PPV}$', color='black')
    ax[0].tick_params(axis='y', labelcolor='black')
    ax[0].set_ylim(0, 1.1)
    ax[0].set_yticks(np.arange(0, 1.1, 0.1))
    ax[1].set_xlabel('Z-score cutoff')
    ax[1].set_ylabel('$\overline{N}_{pairs}$')
    # ax2.set_ylim(0, int(max(data_2[-1]))+10)
    # ax2.set_yticks(range(0, int(max(data_2[-1]))+10, 10))
    ax[1].semilogy()
    plt.xticks(range(int(min(thresh_vals)), int(thresh_vals[-1] + 1)))
    # plt.title(f"{pfam_id} <neff>={int(avg_neff_array[model_n])} (rep{replicate},sub{model_n})")
    ax[0].grid(which="both", alpha=0.3, c="gray")
    ax[1].grid(which="both", alpha=0.3, c="gray")
    for k in range(nmodels):
        # plot_1 = ax[0].errorbar(x, data_1[k], color=colors[k], marker="^", linestyle="dashed",
        #                         yerr=err[k], capsize=6, label=str_neff_l[k])
        plot_1 = ax[0].plot(thresh_vals, data_1[k], color=colors[k], marker="^", linestyle="dashed",
                            label=str_neff_l[k])
        plot_2 = ax[1].plot(thresh_vals, data_2[k], color=colors[k], marker="s")
    ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.40), ncol=5, fancybox=True)
    # ax[0].legend(loc="best")
    # ax[1].legend(loc="best")
    # imgfile = os.path.join(image_dir, f"err_avg_nap_ppv_z_{dcutoff}A.png")
    imgfile = os.path.join(image_dir, f"avg_nap_ppv_z_{dcutoff}A.png")
    plt.savefig(imgfile, dpi=200, format="png", bbox_inches='tight')
    plt.show()


def plot_error(error, nmodels, str_neff_l, thresh_vals, image_dir, dcutoff):
    cmap = plt.get_cmap("nipy_spectral")
    colors = [cmap(i) for i in np.linspace(0, 1, len(aneff1) + 1)]
    plt.figure(3000000 + nmodels, figsize=(4, 4))
    for ii in range(nmodels):
        plt.plot(error[ii, :], c=colors[ii], label=str_neff_l[ii])
    plt.xticks(range(0, len(thresh_vals) + 1, 2), range(int(min(thresh_vals)), int(thresh_vals[-1] + 1)))
    # plt.title(f"{pfam_id} <neff>={int(avg_neff_array[model_n])} (rep{replicate},sub{model_n})")
    plt.grid(which="both", alpha=0.3, c="gray")
    plt.grid(which="both", alpha=0.3, c="gray")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.40), ncol=5, fancybox=True)
    imgfile = os.path.join(image_dir, f"std_ppv_z_{dcutoff}A.png")
    plt.savefig(imgfile, dpi=200, format="png", bbox_inches='tight')
    plt.show()


def pair_ztraj(y, resi, resj, neff_l_str, avg_neff_array, img_dir):
    yavg = y.mean(axis=2)
    yy = yavg.mean(axis=0)
    plt.figure(10, figsize=(11, 4))
    plt.plot(range(len(neff_l_str)), yy, c="black", marker="s")
    plt.title(f"Pair ({resi},{resj}) zscore change")
    plt.xlabel("$<N_{eff}/L>$")
    plt.ylabel("Z-score")
    plt.xticks(range(len(avg_neff_array)), neff_l_str)
    plt.grid(which="both", alpha=0.4, color="gray")
    imgfile = os.path.join(img_dir, f"change_z_pair_{resi}_{resj}.png")
    plt.savefig(imgfile, dpi=200, format="png", bbox_inches='tight')
    plt.show()


DIR_ROOT = "D:\\bootstrap\\"
sys_file = os.path.join(DIR_ROOT, "systems.csv")
df_systems = pd.read_csv(sys_file, header=0)


# pdb = df_systems[df_systems["pfam_id"] == pfam_id].pdb_id[0]


def run_nap(pfam_id):
    DIR_SYS = os.path.join(DIR_ROOT, pfam_id)
    file_filtered_msa = os.path.join(DIR_SYS, f"{pfam_id}_full_filtered_25.fasta")
    L = check_length(file_filtered_msa)

    DIR_RESULTS = os.path.join(DIR_SYS, "results")
    DIR_AVG_RESULTS = os.path.join(DIR_RESULTS, "average_ppv")
    DIR_PDB = os.path.join(DIR_SYS, "PDB")
    DIR_REPLICATES = os.path.join(DIR_SYS, "replicates")
    img_dir = DIR_AVG_RESULTS

    # dca_file =
    imgdir = os.path.join(DIR_AVG_RESULTS, "replicate_ppv_distribution")
    neff_file = os.path.join(DIR_REPLICATES, "neff_array.npy")
    neff_array = np.load(neff_file)
    avg_neff_array = np.mean(neff_array, axis=0)

    nreps = 100
    n_models = avg_neff_array.shape[0]
    d_cutoff = 8
    z_values = np.linspace(-1, 15, 32)
    n_values = np.arange(10, 110, 10)
    ppv = np.zeros((nreps, n_models, len(z_values)))
    total_pairs = np.zeros_like(ppv)
    y = np.zeros_like(ppv)
    for rep in range(nreps):
        for model in range(n_models):
            dir_dca_results = f"{DIR_REPLICATES}\\sub{rep}\\mf\\pc0.2"
            # img_dir = os.path.join(DIR_AVG_RESULTS, "images")
            neff = int(neff_array[rep][model])
            outfile = os.path.join(dir_dca_results, f"{pfam_id}_neff{neff}_pc0.2_all.txt")
            df_dca = pd.read_csv(outfile, delimiter="\t")
            # resi = 24
            # resj = 31
            # y[rep, model] = df_dca.loc[(df_dca["i"] == resi) & (df_dca["j"] == resj)].zscore.values[0]

            # count number of pairs
            tp_count = np.zeros_like(z_values)
            fp_count = np.zeros_like(z_values)
            # total_pairs = np.zeros_like(z_values)
            for j in range(len(z_values)):
                z = z_values[j]
                df_zfilter = df_dca[df_dca["zscore"] >= z].reset_index(drop=True)
                tp_count[j] = len(df_zfilter[df_zfilter['d'] <= d_cutoff])
                fp_count[j] = len(df_zfilter) - tp_count[j]
                total_pairs[rep, model, j] = len(df_zfilter)

            if (total_pairs[rep, model] > 0).all():
                ppv[rep, model] = tp_count / total_pairs[rep, model]
            else:
                ppv[rep, model] = 0

    # pair_ztraj(y, resi, resj, neff_l_str, avg_neff_array, img_dir)
    data_1 = ppv.mean(axis=0)
    data_2 = total_pairs.mean(axis=0)
    error = ppv.std(axis=0)
    return data_1, data_2, avg_neff_array, L, img_dir, n_models, error


pfam_id = "PF00014"
pf2 = "PF00018"
zvals = np.linspace(-1, 15, 32)
n_z = len(zvals)
d1, d2, aneff1, L1, imgdir1, nmodels1, err1 = run_nap(pfam_id)
d12, d22, aneff2, L2, imgdir2, nmodels2, err2 = run_nap(pf2)
neff_l_str1 = [str(f"{kk / L1:.3f}") for idk, kk in enumerate(aneff1)]
neff_l_str2 = [str(f"{kk / L2:.3f}") for idk, kk in enumerate(aneff2)]

# PLOT ERROR
# plot_error(err1, nmodels1, neff_l_str1, zvals, imgdir1, 8)
# plot_error(err2, nmodels1, neff_l_str1, zvals, imgdir2, 8)
# plot_nap_tp(d1, d2, 0, neff_l_str1, zvals, nmodels1, imgdir1, 8)
# plot_nap_tp(d12, d22, 0, neff_l_str2, zvals, nmodels2, imgdir2, 8)

neff_l1 = np.array([kk / L1 for idk, kk in enumerate(aneff1)])
neff_l2 = np.array([kk / L2 for idk, kk in enumerate(aneff2)])

# DEFINE NEFF/L GROUPS
index1_all = np.ravel(np.argwhere(neff_l1 > 0.1))
index1_gt_100 = np.ravel(np.argwhere((neff_l1 > 100)))
index1_lt_100 = np.ravel(np.argwhere((neff_l1 > 10) & (neff_l1 < 100)))
index1_lt_10 = np.ravel(np.argwhere((neff_l1 > 1) & (neff_l1 < 10)))
index1_lt_1 = np.ravel(np.argwhere((neff_l1 > 0.1) & (neff_l1 < 1)))
index2_all = np.ravel(np.argwhere(neff_l2 > 0.1))
index2_gt_100 = np.ravel(np.argwhere((neff_l2 > 100)))
index2_lt_100 = np.ravel(np.argwhere((neff_l2 > 10) & (neff_l2 < 100)))
index2_lt_10 = np.ravel(np.argwhere((neff_l2 > 1) & (neff_l2 < 10)))
index2_lt_1 = np.ravel(np.argwhere((neff_l2 > 0.1) & (neff_l2 < 1)))

# COMBINE FAMILIES
ppv_neff_all = np.concatenate((d1[index1_all], d12[index2_all]))
ppv_neff_100 = np.concatenate((d1[index1_gt_100], d12[index2_gt_100]))
ppv_neff_lt_100 = np.concatenate((d1[index1_lt_100], d12[index2_lt_100]))
ppv_neff_lt_10 = np.concatenate((d1[index1_lt_10], d12[index2_lt_10]))
ppv_neff_lt_1 = np.concatenate((d1[index1_lt_1], d12[index2_lt_1]))

avg_ppv = [0.5, 0.6, 0.7, 0.8, 0.9]
prob_ppv_all = np.zeros((len(avg_ppv), n_z))
prob_ppv_100 = np.zeros_like(prob_ppv_all)
prob_ppv_lt100 = np.zeros((len(avg_ppv), n_z))
prob_ppv_lt_10 = np.zeros((len(avg_ppv), n_z))
prob_ppv_lt1 = np.zeros((len(avg_ppv), n_z))
for apv in range(len(avg_ppv)):
    for i in range(n_z):
        prob_ppv_all[apv, i] = len(ppv_neff_all[:, i][ppv_neff_all[:, i] > avg_ppv[apv]]) / len(ppv_neff_all[:, i])
        prob_ppv_100[apv, i] = len(ppv_neff_100[:, i][ppv_neff_100[:, i] > avg_ppv[apv]]) / len(ppv_neff_100[:, i])
        prob_ppv_lt100[apv, i] = len(ppv_neff_lt_100[:, i][ppv_neff_lt_100[:, i] > avg_ppv[apv]]) / len(
            ppv_neff_lt_100[:, i])
        prob_ppv_lt_10[apv, i] = len(ppv_neff_lt_10[:, i][ppv_neff_lt_10[:, i] > avg_ppv[apv]]) / len(
            ppv_neff_lt_100[:, i])
        prob_ppv_lt1[apv, i] = len(ppv_neff_lt_1[:, i][ppv_neff_lt_1[:, i] > avg_ppv[apv]]) / len(ppv_neff_lt_100[:, i])
# j = 4
for j in range(len(avg_ppv)):
    ls = "solid"
    plt.figure(j+1000, figsize=(4, 4))
    plt.plot(range(len(zvals)), prob_ppv_all[j, :], label=None, c="black")
    plt.scatter(range(len(zvals)), prob_ppv_all[j, :], edgecolors="black", c="black", label="$N_{eff}/L$>0.1")
    plt.plot(range(len(zvals)), prob_ppv_100[j, :], label=None, c="red")
    plt.scatter(range(len(zvals)), prob_ppv_100[j, :], edgecolors="black", c="red", label="$N_{eff}/L$>100")
    plt.plot(range(len(zvals)), prob_ppv_lt100[j, :], label=None, linewidth=2, linestyle=ls)
    plt.scatter(range(len(zvals)), prob_ppv_lt100[j, :], edgecolors="black", label="$10<N_{eff}/L<100$")
    plt.plot(range(len(zvals)), prob_ppv_lt_10[j, :], label=None, linewidth=2, linestyle=ls)
    plt.scatter(range(len(zvals)), prob_ppv_lt_10[j, :], edgecolors="black", label="$1<N_{eff}/L<10$")
    plt.plot(range(len(zvals)), prob_ppv_lt1[j, :], label=None, linewidth=2, linestyle=ls)
    plt.scatter(range(len(zvals)), prob_ppv_lt1[j, :], edgecolors="black", label="$0.1<N_{eff}/L<1$")

    plt.ylim(-0.1, 1.1)
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.xticks(range(0, len(zvals) + 1, 2), range(int(min(zvals)), int(zvals[-1] + 1)))
    plt.xlabel('Z-score cutoff')
    plt.ylabel("$P(\overline{PPV} \geq n \mid Z)$")
    plt.title("$\overline{PPV} \geq$"+f"{avg_ppv[j]}")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.30), ncol=5, fancybox=True)
    plt.grid(which="both", c="gray", alpha=0.2)
    imgfile = os.path.join(DIR_ROOT, f"BPTI_SH3_p_ppv_gt{avg_ppv[j]}_z_8A.png")
    plt.savefig(imgfile, dpi=200, format="png", bbox_inches='tight')
    plt.show()

# NPAIRS
ncorrect_all = np.concatenate((d2[index1_all], d22[index2_all])) * ppv_neff_all
ncorrect_100 = np.concatenate((d2[index1_gt_100], d22[index2_gt_100])) * ppv_neff_100
ncorrect_lt100 = np.concatenate((d2[index1_lt_100], d22[index2_lt_100])) * ppv_neff_lt_100
ncorrect_lt10 = np.concatenate((d2[index1_lt_10], d22[index2_lt_10])) * ppv_neff_lt_10
ncorrect_lt1 = np.concatenate((d2[index1_lt_1], d22[index2_lt_1])) * ppv_neff_lt_1

for ppv_i in range(len(avg_ppv)):
    # ppv_i = 2
    y_all = prob_ppv_all[ppv_i] * ncorrect_all.mean(axis=0)
    y_100 = prob_ppv_100[ppv_i] * ncorrect_100.mean(axis=0)
    y_lt100 = prob_ppv_lt100[ppv_i] * ncorrect_lt100.mean(axis=0)
    y_lt10 = prob_ppv_lt_10[ppv_i] * ncorrect_lt10.mean(axis=0)
    y_lt1 = prob_ppv_lt1[ppv_i] * ncorrect_lt1.mean(axis=0)

    plt.figure(700000+ppv_i, figsize=(4, 4))
    # plt.plot(range(len(zvals)), y_all, label="$N_{eff}/L>0.04$", c='black', marker="o", linewidth=2)
    # plt.plot(range(len(zvals)), y_100, label="$N_{eff}/L>100$", c='red', marker="o", linewidth=2)
    # plt.plot(range(len(zvals)), y_lt100, label="$10<N_{eff}/L<100$", marker="o")
    # plt.plot(range(len(zvals)), y_lt10, label="$1<N_{eff}/L<10$", marker="o")
    # plt.plot(range(len(zvals)), y_lt1, label="$N_{eff}/L<1$", marker="o")
    plt.plot(range(len(zvals)), y_all, label="$N_{eff}/L>0.04$", c='black', linewidth=2)
    plt.plot(range(len(zvals)), y_100, label="$N_{eff}/L>100$", c='red', linewidth=2)
    plt.plot(range(len(zvals)), y_lt100, label="$10<N_{eff}/L<100$")
    plt.plot(range(len(zvals)), y_lt10, label="$1<N_{eff}/L<10$")
    plt.plot(range(len(zvals)), y_lt1, label="$N_{eff}/L<1$")
    # plt.semilogy()
    plt.ylim(0, 80)
    plt.yticks(range(0, 80, 5))
    plt.xticks(range(0, len(zvals) + 1, 2), range(int(min(zvals)), int(zvals[-1] + 1)))
    plt.xlabel("Z-score cutoff")
    plt.ylabel("$\overline{TP}$")
    plt.title("$\overline{PPV} \geq$"+f"{avg_ppv[ppv_i]}")
    plt.grid(which="both", c="gray", alpha=0.2)
    imgfile = os.path.join(DIR_ROOT, f"lines_BPTI_SH3_avgTP_at_ppv{avg_ppv[ppv_i]}_z_8A.png")
    plt.savefig(imgfile, dpi=200, format="png", bbox_inches='tight')
    plt.legend(loc="best")
    plt.show()