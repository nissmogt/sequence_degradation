import matplotlib.pyplot as plt
import os
import numpy as np


def contact_map_single(couplings=None, monomer=None, n=10, x="text", distance_cutoff=6, symmetric=True):
    if len(couplings.columns) == 0 and len(monomer.columns) == 0:
        raise ValueError("Need to specify at least one of couplings or monomer")
    plt.style.use('./styles/contact_map_single.mplstyle')
    plt.figure(num=np.random.randint(10000))
    if len(monomer.columns) > 0:
        plt.scatter("i", "j", data=monomer[monomer["d"] <= distance_cutoff], label="monomer")
    plt.scatter("i", "j", data=couplings[:n], label=f"{x}")
    plt.legend(loc="best")
    plt.xlabel("resi")
    plt.ylabel("resj")
    plt.show()


def compare_contact_map(dca_dataframe, dca_filtered, pdb_dataframe, pdbid, distance_cutoff, n_effective, zcut,
                        Lseq, index_shift, dir_fig, extra_text=None):
    xy_lim = np.abs(Lseq + index_shift)
    top_2n = 2 * Lseq
    plt.figure(figsize=(6, 6))
    plt.scatter('i', 'j', data=pdb_dataframe[pdb_dataframe["d"] <= distance_cutoff],
                color='lightgrey', s=85, label=f"{pdbid} {distance_cutoff}$\AA$")
    plt.scatter('j', 'i', data=pdb_dataframe[pdb_dataframe["d"] <= distance_cutoff],
                color='lightgrey', s=85, label=None)
    plt.scatter('i', 'j', data=dca_dataframe[:top_2n], label=f"top {top_2n} dca",
                linewidths=1.8, color="red", marker="x", s=25)
    plt.scatter('j', 'i', data=dca_filtered, label=f"$z>=${zcut}", color="blue",
                linewidths=1.8, marker="x", s=25)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=4, fancybox=True)
    plt.xlabel("residue i")
    plt.ylabel("residue j")
    plt.xlim(0, Lseq + 5)
    plt.ylim(0, Lseq + 5)
    plt.title(f"<Neff>={n_effective:.2f}, <Neff>/L={n_effective / Lseq:.2f}")
    imgfile = os.path.join(dir_fig, f"{extra_text}_z{zcut}_top{top_2n}_{n_effective:.2f}.png")
    plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_score_distribution(dca_score, n_seqs, dca_dataframe, n_effective, n_res, img_dir, extra_text=None):
    plt.figure(n_seqs, figsize=(5, 5))
    heights, bins, patches = plt.hist(dca_dataframe[dca_score], bins=50)
    plt.semilogy()
    plt.title(f"Neff={n_effective:.2f}, Neff/L={n_effective / n_res:.2f}")
    plt.xlim(-10, 21)
    plt.xticks(np.arange(-10, 24, 3))
    mean = dca_dataframe[dca_score].mean()
    std = dca_dataframe[dca_score].std()
    plt.vlines(mean, 0, 1000, color="black", linestyles="dashed", label=f"mean={mean:.3f}, std={std:.3f}")
    plt.xlabel("Z-score")
    plt.ylabel("Counts")
    plt.legend(loc="upper right")
    img_out = os.path.join(img_dir, f"{extra_text}_{dca_score}_neff{n_effective:.2f}.png")
    plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_dist_distribution(dca_dataframe, n_effective, n_res, img_dir, extra_text=None):
    plt.figure(8927384929, figsize=(5, 5))
    heights, bins, patches = plt.hist(dca_dataframe["d"], bins=50, color="red", edgecolor="black")
    plt.ylim(0, 100)
    plt.title(f"Neff={n_effective:.2f}, Neff/L={n_effective / n_res:.2f}")
    plt.xlabel("intra-chain distance $\AA$")
    plt.ylabel("counts")
    plt.legend(loc="best")
    img_out = os.path.join(img_dir, f"{extra_text}_distance_neff{n_effective:.2f}.png")
    plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_ppv(zfilter, ppv_list, n_effective, pfamid, sequence_len, threshold_val, img_dir, extra_text=None):
    plt.figure(figsize=(5, 5))
    plt.plot(ppv_list)
    plt.scatter(range(len(ppv_list)), ppv_list, edgecolors="black", label=f"z={threshold_val:.1f}")
    plt.ylim(-0.1, 1.1)
    plt.xlabel("rank")
    plt.ylabel("PPV")
    plt.title(f"{pfamid}, Neff:{n_effective:.2f}, Neff/L={n_effective / sequence_len:.2f}")
    plt.grid(which="both", alpha=0.3)
    plt.legend(loc="best")

    if zfilter:
        img_out = os.path.join(img_dir, f"{extra_text}_ppv_neff{n_effective:.2f}_z{threshold_val:.1f}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    else:
        img_out = os.path.join(img_dir, f"{extra_text}_ppv_neff{n_effective:.2f}_top{threshold_val:.1f}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_top_zscore(dca_dataframe, n, n_effective, n_res, img_dir, extra_text=None):
    topn_dca = dca_dataframe["zscore"][:n]
    avg_z = np.mean(topn_dca)
    std_z = np.std(topn_dca)
    plt.figure(873499, figsize=(5, 5))
    plt.plot(topn_dca)
    plt.scatter(range(n), topn_dca, edgecolors="black", label=f"avg:{avg_z:.2f}, std:{std_z:.2f}")
    plt.title(f"Neff={n_effective:.2f}, Neff/L={n_effective / n_res:.2f}")
    plt.xlabel("rank")
    plt.ylabel("z-score")
    plt.ylim(0, 25)
    plt.legend(loc="best")
    img_out = os.path.join(img_dir, f"{extra_text}_top{n}_avgz_neff{n_effective:.2f}.png")
    plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_average_ppv(ppv_array, n_effective_array, sysid, _sys_l, z, _dir_out, norm=1):
    plt.figure(0)
    _r, _n = n_effective_array.shape
    avg_neff_l = np.mean(n_effective_array, axis=0) / _sys_l
    avg_ppv = np.mean(ppv_array, axis=0)
    std_ppv = np.std(ppv_array, axis=0)
    outfile = os.path.join(_dir_out, f"avgppv_std_z{z}.png")
    plt.errorbar(avg_neff_l, avg_ppv, yerr=std_ppv, capsize=6)
    plt.scatter(avg_neff_l, avg_ppv * norm, label=f"z:{z}")
    plt.ylim(-0.1, 1.1)
    plt.legend(loc="best")
    plt.title(f"{sysid}, nreps:{_r}")
    plt.xlabel("average Neff/L")
    plt.ylabel("average ppv")
    plt.semilogx()
    plt.grid(which="both", alpha=0.2)
    plt.savefig(outfile, format="png", dpi=200, bbox_inches='tight')
    plt.close()
    # plt.show()


def multiple_plot_average_ppv(ppv_array, n_effective_array, sysid, _sys_l, z_list, _dir_out, norm=1, extra_text=""):
    _r, _n = n_effective_array.shape
    avg_neff_l = np.mean(n_effective_array, axis=0) / _sys_l
    plt.figure()
    outfile = os.path.join(_dir_out, f"{extra_text}_avgppv_std_z.png")
    for z, z_val in enumerate(z_list):
        avg_ppv = np.mean(ppv_array[z], axis=0)
        std_ppv = np.std(ppv_array[z], axis=0)
        plt.errorbar(avg_neff_l, avg_ppv, yerr=std_ppv, capsize=6)
        plt.scatter(avg_neff_l, avg_ppv * norm, label=f"z:{z_val}")
    plt.ylim(-0.1, 1.1)
    plt.legend(loc="best")
    plt.title(f"{sysid}, nreps:{_r}")
    plt.xlabel("average Neff/L")
    plt.ylabel("average ppv")
    plt.semilogx()
    plt.grid(which="both", alpha=0.2)
    plt.savefig(outfile, format="png", dpi=200, bbox_inches='tight')
    plt.close()
    # plt.show()


def plot_avg_zscore(z_array, n_effective, pdbid, n_res, img_dir, extra_text=""):
    from analysis.zscore import calculate_average_zscore
    average_z, std_z = calculate_average_zscore(z_array)
    plt.figure(7399, figsize=(5, 5))
    _m, _n = average_z.shape
    neff_l = n_effective / n_res
    for i in range(_n):
        plt.plot(average_z[i, :])
        plt.scatter(range(_n), average_z[i, :], edgecolors="black", label=f"neff/L:{neff_l[i]:.2f}")
    plt.title(f"{pdbid}")
    plt.xlabel("di rank")
    plt.ylabel("average z-score")
    plt.ylim(0, 25)
    plt.legend(loc="best")
    img_out = os.path.join(img_dir, f"{extra_text}_avgz_top{_n}.png")
    # plt.show()
    plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_avg_dist(dist_array, cutoff, n_effective, pdbid, n_res, img_dir, extra_text=""):
    average_dist = np.mean(dist_array, axis=0)
    plt.figure(799, figsize=(6, 8))
    _m, _n = average_dist.shape
    neff_l = n_effective / n_res
    for i in range(_n):
        plt.plot(range(1, _n+1), average_dist[i, :])
        plt.scatter(range(1, _n+1), average_dist[i, :], edgecolors="black", label=f"neff/L:{neff_l[i]:.2f}")
    plt.fill_between(range(1, _n+1), cutoff, alpha=0.2)
    plt.hlines(cutoff, 1, _n, colors="black", linestyles="dashed", label=f"{cutoff}Å")
    plt.title(f"{pdbid}")
    plt.xlabel("di rank")
    plt.ylabel("average distance (Å)")
    plt.ylim(0, 30)
    plt.legend(loc="best")
    img_out = os.path.join(img_dir, f"{extra_text}_avg_dist_top{_n}.png")
    # plt.show()
    plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()
