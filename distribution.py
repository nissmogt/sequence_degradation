import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# f = "contact_maps\\atom_distance_matrix_7lvs.txt"
# f = "contact_maps\\atom_distance_matrix_1R8U.txt"


def cut(df_inter, range_1):
    # Cut the fusion peptide into individual segments
    seg_1 = 15
    seg_2 = seg_1 + 4
    seg_3 = seg_2 + 31

    # check the distances between cited and taz
    dist_cited = []
    for i, elem in (df_inter.iterrows()):
        for k, idx in enumerate(range_1[seg_2:seg_3]):
            # in this case, j is the residue from desired chain (i.e., B)
            if (elem.j == idx):
                dist_cited.append([elem.i, elem.j, elem.si, elem.sj, elem.d])
    # TODO: write routine for other pdbs and compare distributions
    return np.array(dist_cited)


def peptide_cut(df_inter, range_1):
    # Cut the fusion peptide into individual segments
    seg_cited = 22
    seg_linker = seg_cited + 4
    seg_hif = seg_linker + 30

    # check the distances between cited and taz
    dist_cited = []
    for i, elem in (df_inter.iterrows()):
        for k, idx in enumerate(range_1[:seg_cited]):
            # for k, idx in enumerate(range_1[:seg_cited]):
            # in this case, j is the residue from desired chain (i.e., B)
            if elem.j == idx:
                dist_cited.append([elem.i, elem.j, elem.si, elem.sj, elem.d])
    # TODO: write routine for other pdbs and compare distributions
    return np.array(dist_cited)


def function_(f, d_threshold, flag):
    """

    :param flag:
    :param d_threshold:
    :param f:
    :return:
    """
    # d_threshold = 8  # contact definition
    df = pd.read_csv(f, delimiter="\t", header=0)
    df_mon = df[df.chain_1 == df.chain_2]
    inter = df[df.chain_1 != df.chain_2]
    # select domain
    cf = max(df_mon[df_mon.chain_1 == "B"].j)  # get the end of entire protein
    ci = min(df_mon[df_mon.chain_1 == "B"].i)  # get the end of chain A
    range_n = np.arange(ci, cf)
    inter = inter[inter.d < d_threshold]  # contacts are defined by d_threshold
    if flag:
        return cut(inter, range_n)
    else:
        return peptide_cut(inter, range_n)


d = 12
f1 = "contact_maps\\atom_distance_matrix_7lvs.txt"
# f3 = "contact_maps\\atom_distance_matrix_1L8C.txt"
f3 = "contact_maps\\atom_distance_matrix_1R8U.txt"
pdbid = f3.split("\\")[1].strip(".txt")[-4:]
outfile = "7lvs_{}_difference.txt".format(pdbid.upper())
trimer_pairs = function_(f1, d, False)
# dimer_pairs = function_(f3, d, True)
dimer_pairs = function_(f3, d, False)

# reindex 1R8U i,j bc their chains are not in 7LVS order
# for k, kdx in enumerate(trimer_pairs):
# trimer_pairs[k] = [kdx[0], kdx[1] - 3, kdx[2], kdx[3], kdx[-1]]
# For 1R8U
for i, idx in enumerate(dimer_pairs):
    dimer_pairs[i] = [idx[0] - 1, idx[1] - 3, idx[2], idx[3], idx[-1]]
# For 1L8C
# for i, idx in enumerate(dimer_pairs):
#     dimer_pairs[i] = [idx[0] + 4, idx[1] + 9, idx[2], idx[3], idx[-1]]

om = []
gilgamesh = []
jchrist = []
moham = []
neo = []
# categorize delta_td points into 3 classes for easy distinction
for i, idx in enumerate(trimer_pairs):
    for j, jdx in enumerate(dimer_pairs):
        if idx[0] == jdx[0] and idx[1] == jdx[1]:
            rij_t = idx[-1]
            rij_d = jdx[-1]
            delta_td = rij_t - rij_d
            om.append([idx[0], idx[1], idx[2], idx[3], jdx[2], jdx[3], delta_td, rij_t, rij_d])
            if rij_t <= 5.4:
                gilgamesh.append([idx[0], idx[1], idx[2], idx[3], jdx[2], jdx[3], delta_td, rij_t, rij_d])
            if rij_d <= 5.4:
                jchrist.append([idx[0], idx[1], idx[2], idx[3], jdx[2], jdx[3], delta_td, rij_t, rij_d])
            if rij_t > 5.4:
                moham.append([idx[0], idx[1], idx[2], idx[3], jdx[2], jdx[3], delta_td, rij_t, rij_d])
            if rij_d > 5.4:
                neo.append([idx[0], idx[1], idx[2], idx[3], jdx[2], jdx[3], delta_td, rij_t, rij_d])
om = np.array(om)
header = ["i", "j", "t_i", "t_j", "d_i", "d_j", "delta_td", "rij_t", "rij_d"]
df_om = pd.DataFrame(om, columns=header)
df_g = pd.DataFrame(gilgamesh, columns=header)
df_j = pd.DataFrame(jchrist, columns=header)
x_trimer = df_om[df_om.delta_td <= -5]
x_dimer = df_om[df_om.delta_td >= 5]
# df_om.to_csv("7lvs_1L8C_difference.csv", header=["i", "j", "t_i", "t_j", "d_i", "d_j", "delta_td", "rij_t", "rij_d"])
df_om['i'] = df_om['i'].astype(int)
df_om['j'] = df_om['j'].astype(int)
df_om['t_i'] = df_om['t_i'].astype(int)
df_om['t_j'] = df_om['t_j'].astype(int)
df_om['d_i'] = df_om['d_i'].astype(int)
df_om['d_j'] = df_om['d_j'].astype(int)
df_om.to_csv(outfile, sep='\t', index=False, header=header, float_format="%.4f")

# plt.figure(0)
# plt.hist(trimer_pairs[:, -1], label="ternary structure")
# plt.hist(dimer_pairs[:, -1], label="dimer", alpha=0.6)
# plt.ylabel("counts")
# plt.xlabel("distance (A)")
# plt.legend(loc="best")
# plt.show()

# Z-score
# mean = np.mean(df_om.delta_td)
# std = np.std(df_om.delta_td)
# zscore = np.zeros(len(df_om))
# for i, idx in enumerate(df_om.delta_td):
#     zscore[i] = (idx - mean) / std

# plt.figure(10)
# plt.scatter(range(len(zscore)), zscore)
# plt.ylabel("Z-score")
# plt.xlabel("array index")
# plt.show()

# plt.figure(1)
# plt.scatter(range(len(df_om)), df_om.delta_td)
# plt.hlines(4, xmin=0, xmax=len(df_om), color="black", label="")
# plt.hlines(-4, xmin=0, xmax=len(df_om), color="black", label="")
# plt.hlines(2.92, xmin=0, xmax=len(df_om), color="black", label="$\sigma_{t}$+$\sigma_{t}$")
# plt.hlines(-2.92, xmin=0, xmax=len(df_om), color="black")
# plt.hlines(2.13, xmin=0, xmax=len(df_om), color="red", label="$\sigma_{t}$*$\sigma_{t}$")
# plt.hlines(-2.13, xmin=0, xmax=len(df_om), color="red")
# plt.xlabel("sample points")
# plt.ylabel("$\Delta_{td}$")
# plt.legend(loc="best")
# plt.show()
