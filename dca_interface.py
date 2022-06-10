import os.path
import pandas as pd
import matplotlib.pyplot as plt
from rank_hamming import rank_hamming
import numpy as np
from Bio import pairwise2
from mapping_functions import map_indices, apply_map
from plot_contact_map import plot_cm
from zscore import zscore_calc
from output_pipeline import (pipeline_pdb_distance_matrix, pipeline_interface)
from pdb_index_converter import align
from universal_reference_distribution import ref_dist


def tp(df, n, dcut):
    count = 0
    ncount = 0
    df = df.reset_index(drop=True)
    _df = df[:n]
    tp = np.zeros(len(_df))
    fp = np.zeros(len(_df))
    for ii in range(len(_df)):
        if _df.d[ii] <= dcut:
            count += 1
        else:
            ncount += 1
        tp[ii] = count
        fp[ii] = ncount
    plt.figure(np.random.randint(11, 20))
    plt.plot(range(len(tp)), tp, label="True Positive, 6A", marker="o")
    plt.plot(range(len(fp)), fp, label="False Positive", marker="o")
    plt.plot(range(len(fp)), range((len(fp))), linestyle="dashed", label="x=y", color="black")
    plt.xlabel("rank index")
    plt.ylabel("count")
    plt.legend(loc="best")
    plt.show()


def contact_map(df1, df2, df3, N, title):
    plt.figure(np.random.randint(10))
    plt.scatter("i", "j", data=df1, color="gray", label="Native Monomer")
    plt.scatter("i", "j", data=df2, color="black", label="Native Dimer", marker="s")
    plt.scatter("i", "j", data=df3[:N], color="red", marker="s", label="Top {} DCA pred.".format(N),
                edgecolors='black')
    plt.legend(loc="best")
    plt.xlabel("residue i")
    plt.ylabel("residue j")
    plt.title(title)
    plt.show()


def histogram(x, l):
    # plt.figure(np.random.randint(10, 20))
    plt.hist(x, bins=100, label=l)
    plt.semilogy()
    plt.xlabel("CN-score")
    plt.ylabel("Counts")
    plt.legend(loc="best")
    plt.show()


def make_lists(df_dca, difference_list, thresh):
    # make a list of unique dimer, trimer, monomer, and unassigned contacts
    trimer_thresh = thresh[0]
    dimer_thresh = thresh[1]
    trimer_list = []
    dimer_list = []
    both_list = []
    unassigned_list = []
    # union of z-score-filtered pairs with unique list
    for i, elem in df_dca.iterrows():
        u = elem.si
        v = elem.sj
        cn = elem.fn_apc
        zscore = elem.zscore_calculation
        for j, elem2 in difference_list.iterrows():
            k = int(elem2.t_i)
            m = int(elem2.t_j)
            if u == elem2.d_i and v == elem2.d_j:
                if elem2.delta_td < trimer_thresh:
                    trimer_list.append([u, v, k, m, cn, elem2.rij_t, elem2.rij_d,
                                        elem2.delta_td, zscore])
                if elem2.delta_td > dimer_thresh:
                    dimer_list.append([u, v, k, m, cn, elem2.rij_t, elem2.rij_d,
                                       elem2.delta_td, zscore])
                else:
                    both_list.append([u, v, k, m, cn, elem2.rij_t, elem2.rij_d,
                                      elem2.delta_td, zscore])
            else:
                unassigned_list.append([u, v, k, m, cn, elem2.rij_t, elem2.rij_d,
                                        elem2.delta_td, zscore])

    h = ["i", "j", "t_i", "t_j", "cn", "rij_t", "rij_d", "delta_td", "zscore"]
    df_trimer = pd.DataFrame(trimer_list, columns=h)
    df_dimer = pd.DataFrame(dimer_list, columns=h)
    df_both = pd.DataFrame(both_list, columns=h)
    df_u = pd.DataFrame(unassigned_list, columns=h)
    return [df_trimer, df_dimer, df_both, df_u]


def mapping(map_dca2pdb_dict, ranked_dca, pdb_df_list):
    mapped_dca_array = apply_map(ranked_dca.to_numpy(), map_dca2pdb_dict)
    # df_dca_mapped = pd.DataFrame(mapped_dca_array, columns=['i', 'j', 'fn_apc', 'fn', 'ui', 'uj'])
    df_dca_mapped = pd.DataFrame(mapped_dca_array, columns=['i', 'j', 'di', 'di_apc', 'mi', 'ui', 'uj'])
    df_dca_mapped['i'] = df_dca_mapped['i'].astype(int)
    df_dca_mapped['j'] = df_dca_mapped['j'].astype(int)
    df_dca_mapped['ui'] = df_dca_mapped['ui'].astype(int)
    df_dca_mapped['uj'] = df_dca_mapped['uj'].astype(int)
    df_dca_mapped_dist = df_dca_mapped.merge(pdb_df_list[0], how='inner', on=['i', 'j'])
    return df_dca_mapped_dist


def index_map(outfile, msa_seq, pdb_seq):
    # need to penalize for opening and adding gaps otherwise mapping is off (s param {-.5,-.1})
    # Index mapping
    reference_alignment = align(msa_seq, pdb_seq)
    # alignments = pairwise2.align.globalxs(pdb_seq, msa_seq, -.5, -.1)
    # print(pairwise2.format_alignment(*alignments[0]))
    # m = map_indices(alignments[0][0], 1, 0, alignments[0][1], 1, 0)
    reference_alignment = reference_alignment.rename(columns={"i": "pdb_i", "A_i": "pdb_res",
                                                              "j": "dca_i", "A_j": "dca_res"})
    np.savetxt(outfile, reference_alignment, header="pdb_i\tpdb_res\tdca_i\tdca_res", fmt="%s\t%s\t%s\t%s",
               comments='')
    print("(map_dca2pdb)\tWrote {}".format(outfile))
    map_pdb_dca_list = reference_alignment.dropna()
    # convert mapping list to dictionary for easier referencing
    map_dca2pdb_dict = dict(zip(map_pdb_dca_list["dca_i"], map_pdb_dca_list["pdb_i"]))
    return map_dca2pdb_dict


def z_score(df, zoutfile, calc, ref, output_ref_dist=False):
    """

    :param df: dca dataframe
    :param zoutfile: outfile name
    :param calc: if true, calculates zscore
    :param ref: if true, makes a reference distribution
    :return:
    """
    if ref:
        reference_distribution = ref_dist()
    else:
        reference_distribution = np.load("zref.npy")

    if calc:
        # zheader = ["i", "j", "fn_apc", "fn", "ui", "uj", "d", "si", "sj", "chain_1", "chain_2", "resnames",
        #            "atom_id", "zscore"]
        zheader = ["i", "j", "di", "diapc", "mi", "d", "si", "sj", "chain_1", "chain_2", "resnames",
                   "atom_id", "zscore"]
        df_zscore = zscore_calc(df, ref=reference_distribution, which='all', score='di')
        df_zscore.to_csv(zoutfile, sep='\t', index=False, header=zheader, float_format='%.4f')

    else:
        df_zscore = pd.read_csv(zoutfile, delimiter="\t")

    if output_ref_dist:
        return df_zscore, reference_distribution
    else:
        return df_zscore


def dca(raw_dca, score, template_alignment, dca_out, df_header, reference_outfile, pdb_df, align=False):
    """
    Load DCA output, sort by score, map to pdb index, and return interchain DCA pairs
    :param df_header: header to use for df
    :param raw_dca: raw dca output
    :param template_alignment: File with pdb and msa seq to align
    :param dca_out: outfile for mapped and ranked dca
    :param align: boolean; if True, align using template_alignment
    :return:
    """
    raw = pd.read_csv(raw_dca, delimiter=',', names=['i', 'j', 'di', 'di_apc', 'mi'])
    # raw = pd.read_csv(raw_dca, delimiter=',', names=['i', 'j', 'fn_apc', 'fn'])
    dca_ranked = rank_hamming(raw, score, 5)

    # Sequences
    with open(template_alignment, 'r') as fp:
        pheader = fp.readline().rstrip()
        pdb_seq = fp.readline().rstrip()
        mheader = fp.readline().rstrip()
        msa_seq = fp.readline().rstrip()
    assert (len(msa_seq) > 1 and len(pdb_seq) > 1)

    if align:
        # Map dca to pdb index
        map_dca2pdb_dict = index_map(reference_outfile, msa_seq, pdb_seq)
        df_dca_mapped_dist = mapping(map_dca2pdb_dict, dca_ranked, pdb_df)
        df_dca_mapped_dist.to_csv("{}_all_dca.txt".format(dca_out[:4]), sep='\t', index=False, header=df_header,
                                  float_format='%.4f')
        # df_dca_mapped_dist.to_csv("{}_all_dca.txt".format(dca_out[:4]), sep='\t', index=False, header=df_header, float_format='%.4f')
        df_interface = df_dca_mapped_dist[df_dca_mapped_dist["chain_1"] != df_dca_mapped_dist["chain_2"]]
        df_interface.to_csv(dca_out, sep='\t', index=False, header=df_header, float_format='%.4f')
    else:
        df_interface = pd.read_csv(dca_out, delimiter="\t")

    return df_dca_mapped_dist, df_interface


def generate_filenames(id_pdb, raw_dca):
    """
    Helper function that generates filenames needed for other functions.
    :param id_pdb:
    :param raw_dca:
    :return:
    """
    grain = 'aa'
    if grain == 'aa':
        # header = ["i", "j", "fn_apc", "fn", "ui", "uj", "d", "si", "sj", "chain_1", "chain_2", "resnames", "atom_id"]
        header = ["i", "j", "di", "di_apc", "mi", "ui", "uj", "d", "si", "sj", "chain_1", "chain_2", "resnames",
                  "atom_id"]
    else:
        header = ["i", "j", "fn_apc", "fn", "ui", "uj", "d", "si", "sj", "chain_1", "chain_2", "resnames"]

    template_align = "align_{}.fasta".format(id_pdb)
    align_reference = "ref_map_{}.txt".format(raw_dca.strip(".fas.txt"))
    dca_out = "{}_{}_inter_mapped_{}_dist.txt".format(id_pdb, raw_dca.strip(".fas.txt"), grain)
    trimer_difference = "7lvs_{}_difference.txt".format(id_pdb.upper())
    # zscore_outfile = "zscore_" + dca_out
    zscore_outfile = "zscore_" + "{}_all_dca.txt".format(id_pdb)
    return [template_align, align_reference, dca_out, trimer_difference, zscore_outfile, header]


def gg(pdb_in, raw_dca_in, zcut):
    # Params: and Dimer/Trimer difference threshold
    t_thresh = -3
    d_thresh = 3
    #
    pdbid = os.path.basename(pdb_in).strip(".pdb")
    template_alignment_file, reference_outfile, dca_fileout, diff_list, zout, header = generate_filenames(pdbid,
                                                                                                          raw_dca_in)

    # Load/Calculate PDB distance matrix
    pdb_df_list, chain_lengths = pipeline_pdb_distance_matrix(pdb_in, cutoff=12, heavy_atom=True, read=True)
    mon = pdb_df_list[1]
    inter = pdb_df_list[2]

    score = "di_apc"
    # DCA interface
    df_all, df_inter = dca(raw_dca_in, score, template_alignment_file, dca_fileout, header, reference_outfile, pdb_df_list,
                           align=True)
    # plt.figure(1000)
    # plt.scatter("fn_apc", "d", data=df_all, label="All DCA predictions")
    # plt.scatter("fn_apc", "d", data=df_inter, label="DCA Inter-chain predictions", marker="x", color="red")
    # plt.xlabel("FN_apc")
    # plt.ylabel("Distance (${\AA}$)")
    # plt.xlim(xmax=2.5, xmin=-0.2)
    # plt.ylim(ymax=70, ymin=0)
    plt.figure(2000)
    plt.hist("fn_apc", data=df_all, bins=50, edgecolor='black', label="fn-apc")
    plt.hist("fn", data=df_all, bins=50, edgecolor='black', alpha=0.6,label="fn")
    plt.semilogy()
    plt.xlabel("score")
    plt.ylabel("counts")
    plt.legend(loc="best")
    plt.show()

    # calculate Zscore
    # df_z, rd_ = z_score(df_all, zout, calc=True, ref=True, output_ref_dist=True)
    # df_z = z_score(df_inter, zout, calc=True, ref=True)
    # z_filter = df_z[df_z.zscore >= zcut]

    # import unique lists
    # df_list = pd.read_csv(diff_list, header=0, delimiter="\t")
    # trimer, dimer, both, and unassigned lists
    # tdbu = make_lists(z_filter, df_list, [t_thresh, d_thresh])
    contact_map(mon, inter, df_all, 10, pdbid)
    contact_map(mon, inter, df_inter, 10, pdbid)
    # histogram(rd_, l="Zscore reference from paper")
    # tp(df_all, n=10, dcut=6)
    # tp(df_inter, n=10, dcut=6)
    print("random")
    return tdbu


# File inputs
pdb_in_cited = "structures\\1R8U.pdb"
pdb_in_hif = "structures\\1L8C.pdb"
raw_dca_in_cited = "DI_TAZ1_CITED2_pc0.2.di"
# raw_dca_in_cited = "DI_PC-0.80.txt"
# raw_dca_in_hif = "FN_APC_TAZ1_HIF1A_hmmer.fas.txt"

# Z-score cutoff
z_threshold = 5
# Run Program and return unique lists according to a zscore value
# cited = gg(pdb_in_cited, raw_dca_in_cited, z_threshold)
# hif = gg(pdb_in_hif, raw_dca_in_hif, z_threshold)

# Output to lists to csv file
# nm = ["cited", "hif"]
# nm = ["cited"]
# out = []
# for q, label in enumerate(nm):
#     out.append(
#         ["{}_unique_trimer_{}.csv".format(label, z_threshold), "{}_unique_dimer_{}.csv".format(label, z_threshold),
#          "{}_both_{}.csv".format(label, z_threshold), "{}_unassigned_{}.csv".format(label, z_threshold)]
#     )

# l = len(cited)
# assert l == len(hif)
# TODO: Eventually implement output into make_lists()
# for i in range(l):
#     cited[i].to_csv(out[0][i], index=False, float_format='%.4f')
    # hif[i].to_csv(out[1][i], index=False, float_format='%.4f')

# n = 10000
# mon = mon[mon["chain_1"] == 'A']
# plt.scatter('i', 'j', data=mon[mon["d"] <= 8], color='gray')
# plt.scatter('i', 'j', data=inter[inter["d"] <= 12], color='tan')
# plt.scatter('i', 'j', data=df_inter[:n], alpha=1, color='red')
# plt.xlabel("distance")
# plt.ylabel("fn_apc")
# plt.show()
