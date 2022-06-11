import matplotlib.pyplot as plt
from draw_contacts import draw_publish_dca
import os
import seaborn as sns
from collections import defaultdict
import pandas as pd
import numpy as np
from distance_pdb import distance_matrix
from rank_hamming import rank_hamming
from msa_functions import msa_object, effective_seq_calc
import matlab.engine
from gapFilter import gap_filter
from msa_functions import generate_replicates

plt.style.use('paper.mplstyle')

# GLOBAL PARAMETERS
# DIR_ROOT = "D:\\bootstrap\\"
# DIR_DCA_CODE = "C:\\Users\\kmehr\\MATLAB Drive\\Research\\dca_package"
DIR_ROOT = "bootstrap"
DIR_DCA_CODE = "Users\\euler\\MATLAB-Drive\\Research\\Scripts\\mf_plm_reweight"


def make_dir_struct(pfam_id):
    # DIRECTORY STRUCTURE INITIALIZATION
    dir_sys = os.path.join(DIR_ROOT, pfam_id)
    dir_results = os.path.join(dir_sys, "results")
    dir_avg_results = os.path.join(dir_results, "average_ppv")
    dir_pdb = os.path.join(dir_sys, "PDB")
    dir_dist_mat = os.path.join(dir_pdb, "distance_matrix")
    dir_replicates = os.path.join(dir_sys, "replicates")
    dir_nseq_list = os.path.join(dir_replicates, "length_list.txt")

    # MAKE DIRECTORY IF DOESNT EXIST
    list_dir = [dir_sys, dir_results, dir_avg_results, dir_pdb, dir_dist_mat, dir_replicates, dir_nseq_list]
    for entry in list_dir:
        if not os.path.exists(entry):
            os.makedirs(entry)

    return list_dir


def collapse_chain(pdb_dataframe, n_amino):
    # n_amino = 375
    pdb_dataframe.loc[:, ["i", "j"]] = pdb_dataframe.loc[:, ["i", "j"]] % n_amino
    pdb_dataframe["i"].replace(0, n_amino, inplace=True)
    pdb_dataframe["j"].replace(0, n_amino, inplace=True)
    return pdb_dataframe[pdb_dataframe["i"] < pdb_dataframe["j"]]


def into_matrix(dca_dataframe):
    dca_matrix = np.zeros((dca_dataframe["j"].max() + 1, dca_dataframe["j"].max() + 1))
    xx = dca_dataframe.loc[:, ["i", "j", "zscore"]].to_numpy()
    for idx, val in enumerate(xx):
        dca_matrix[int(val[0]), int(val[1])] = val[2]
    return dca_matrix


def reverse_odd_list(olist):
    half_length = int(np.ceil(len(olist) / 2))
    for i in range(half_length - 1):
        tmp = olist[i]
        olist[i] = olist[-(i + 1)]
        olist[-(i + 1)] = tmp
    return olist


def check_length(filein):
    from Bio import AlignIO
    align = AlignIO.read(open(filein), "fasta")
    return align.get_alignment_length()


def pipeline_pdb(PDB_id, dir_pdb):
    # PDB distance matrix
    # PDB_id = "5pti"
    # PDB_id = "1or7"
    dmatrix_file = os.path.join(dir_pdb, f"atom_distance_matrix_{PDB_id}.txt")
    if os.path.exists(dmatrix_file):
        pdb_dataframe = pd.read_csv(dmatrix_file, header=0, delimiter="\t")
    else:
        # uncomment if not reading
        pdb_path = os.path.join(dir_pdb, f"{PDB_id}.pdb")
        pdb_dataframe, chain_len = distance_matrix(pdb_path, heavy_atom=True)

    # min_i = pdb_dataframe[pdb_dataframe.si > 0]["i"].min()
    # pdb_dataframe["i"] = pdb_dataframe["i"] - (min_i - 1)
    # pdb_dataframe["j"] = pdb_dataframe["j"] - (min_i - 1)
    # out_df = pdb_dataframe[pdb_dataframe.si > 0]
    # return out_df.reset_index(drop=True)
    return pdb_dataframe


def pipeline_inference(pfam_id, matlab_input, model_length, pseudocount, theta):
    output_dir = os.path.join(os.path.dirname(matlab_input), f"mf\\pc{pseudocount:.1f}\\")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_dca_output = os.path.join(f"{output_dir}", f"DI_{pfam_id}_n{model_length}.txt")
    file_matrix_output = os.path.join(f"{output_dir}", f"matrix_{pfam_id}_n{model_length}.mat")
    if not os.path.exists(file_dca_output):
        eng = matlab.engine.start_matlab()
        eng.addpath(DIR_DCA_CODE, nargout=0)
        eng.addpath(output_dir, nargout=0)
        neffective = eng.dca_h_J_Full_v4(matlab_input, file_dca_output, file_matrix_output, pseudocount,
                                         theta, nargout=1)
        eng.quit()
        return neffective
    else:
        return 0


def neff_calculation(matlab_input, theta):
    eng = matlab.engine.start_matlab()
    code_dir = "C:\\Users\\kmehr\\MATLAB Drive\\Research\\Scripts\\plm-code"
    function_dir = "C:\\Users\\kmehr\\MATLAB Drive\\Research\\Scripts\\plm-code\\functions"
    eng.addpath(code_dir, nargout=0)
    eng.addpath(function_dir, nargout=0)
    neffective = eng.neff_calc(matlab_input, theta, 6, nargout=1)
    eng.quit()
    return neffective


def run_neff_calculation(pfam_id, len_list, nreplicates, outdir, theta, passthrough=False):
    output = os.path.join(outdir, "neff_array.npy")
    if passthrough:
        # Load neffective array
        return np.load(output)
    else:
        # Calculate neff for every replicate and degraded model
        n_effective_array = np.zeros((nreplicates, len(len_list)))
        for rep in range(nreplicates):
            msa_rep_dir = os.path.join(outdir, f"sub{rep}")
            for model in range(len(len_list)):
                model_length = len_list[model]
                msa_input = os.path.join(msa_rep_dir, f"{pfam_id}_n{model_length}_sub{rep}.fasta")
                print(f"PFAM: {pfam_id} REP: {rep} N{model}: {model_length}")
                n_effective_array[rep][model] = neff_calculation(msa_input, theta)

        # Save Neff to file
        np.save(output, n_effective_array)
        return n_effective_array


def run_replicates(pfam_id, len_list, nreplicates, outdir, pseudocount, theta, passthrough=False):
    output = os.path.join(outdir, "neff_array.npy")
    if passthrough:
        # Load neffective array
        return np.load(output)
    else:
        # Run DCA for every replicate and degraded model
        n_effective_array = np.zeros((nreplicates, len(len_list)))
        for rep in range(nreplicates):
            msa_rep_dir = os.path.join(outdir, f"sub{rep}")
            for model in range(len(len_list)):
                model_length = len_list[model]
                if model_length > 0:
                    msa_input = os.path.join(msa_rep_dir, f"{pfam_id}_n{model_length}_sub{rep}.fasta")
                    print(f"PFAM: {pfam_id} REP: {rep} N{model}: {model_length}")
                    print(f"{msa_input}")
                    n_effective_array[rep][model] = pipeline_inference(pfam_id, msa_input, model_length, pseudocount,
                                                                       theta)

        # Save Neff to file
        np.save(output, n_effective_array)


def zscore_calculation(dca_dataframe, dca_score):
    # calc zscore
    # reference = np.load(f"dimers_{dca_score.upper()}_scores.npy")
    reference = np.load(f"monomer_{dca_score.upper()}_scores.npy")
    scores = dca_dataframe[dca_score].to_numpy()
    new_reference = np.concatenate((reference, scores))
    mean = np.mean(new_reference)
    std = np.std(new_reference)
    z = np.zeros_like(scores)
    for k in range(len(scores)):
        z[k] = (scores[k] - mean) / std
    dca_dataframe.insert(3, "zscore", z)
    return dca_dataframe


def pipeline_process_results(pfam_id, pdb_id, dca_in, dir_system, pc, pdb_dataframe, score_type, n_effective,
                             theta_flag=False, load_processed=False):
    from map_v2 import map_dca
    from map_dca import mapit
    """

    :param theta_flag:
    :param load_processed:
    :param dca_in: Str
    :param pdb_dataframe: DF
    :param score_type: Str
    :param n_effective: Int
    :return:
    """
    df_header = ["i", "j", "di", "zscore", "diapc", "mi", "d", "si", "sj", "chain_1", "chain_2", "resnames", "atom_id"]
    dir_dca_in = os.path.dirname(dca_in)
    if theta_flag:
        outfile = os.path.join(dir_dca_in, f"{pfam_id}_t{n_effective:.2f}_pc{pc}_all.txt")
    else:
        outfile = os.path.join(dir_dca_in, f"{pfam_id}_neff{n_effective}_pc{pc}_all.txt")
    # os.remove(outfile)
    print(outfile)
    if load_processed:
        df_dca_mapped_dist = pd.read_csv(outfile, delimiter="\t")
        index_shift = df_dca_mapped_dist.i[0] - df_dca_mapped_dist.si[0]
    else:
        # LOAD RAW DCA FILE
        df = pd.read_csv(dca_in, delimiter=',', header=0)

        # RANK AND SEQ SEPARATION
        df_rank = rank_hamming(df, score_type, 5)

        # df_mapped = map_dca("2zwh_scan.txt", df_rank)

        # MAPPING TO PDB INDEX
        if pdb_id == "1vhv":
            hmmscan_file = os.path.join(dir_system, f"{pdb_id.lower()}_scan.txt")
            df_mapped = mapit(df_rank, pdb_id, dir_system, align=hmmscan_file)

            # ZSCORE
            df_z = zscore_calculation(df_mapped, score_type)
            index_shift = 0
        else:
            # ZSCORE
            df_z = zscore_calculation(df_rank, score_type)
            mapping_file = os.path.join("C:\\Users\\kmehr\\PycharmProjects\\idp_trimer", "pdb_pfam_mapping.txt")
            df_mapping = pd.read_csv(mapping_file, delimiter="\t", header=0, comment="#")

            pdb_ = df_mapping[df_mapping.PDB == pdb_id]
            pdb_chain = pdb_dataframe.chain_1[0]
            pfam_entry = pdb_[pdb_.PFAM_ACCESSION == pfam_id]
            chain_filter = pfam_entry[pfam_entry.CHAIN == pdb_chain]
            pfam_start = int(chain_filter.AUTH_PDBRES_START)
            pdb_entry = int(chain_filter.PDB_START)

            if pdb_id == "1oap":
                index_shift = pdb_entry - 3
            else:
                index_shift = pdb_entry - 1
            df_z["i"] = df_z["i"] + index_shift
            df_z["j"] = df_z["j"] + index_shift

        df_out = df_z
        df_dca_mapped_dist = df_out.merge(pdb_dataframe, how='inner', on=['i', 'j'])
        df_dca_mapped_dist.to_csv(outfile, sep='\t', index=False, header=df_header, float_format='%.4f')

    return df_dca_mapped_dist, index_shift


def calculate_ppv(dataframe_in, threshold):
    dca_dataframe = dataframe_in.reset_index(drop=True)
    tp_count = np.zeros(len(dca_dataframe))
    fp_count = np.zeros_like(tp_count)
    positive_predictive_value = np.zeros_like(tp_count)
    tcount = 0
    fcount = 0

    for j in range(0, len(dca_dataframe)):
        if dca_dataframe["d"][j] <= threshold:
            tcount += 1
            tp_count[j] = tcount
            fp_count[j] = fcount
        else:
            fcount += 1
            tp_count[j] = tcount
            fp_count[j] = fcount
        positive_predictive_value[j] = tp_count[j] / (fp_count[j] + tp_count[j])

    # print(positive_predictive_value)
    # return tp_count
    return positive_predictive_value, tp_count + fp_count


def figure_1(dca_dataframe, dca_filtered, pdb_dataframe, pdbid, distance_cutoff, average_neff, zcut,
             Lseq, index_shift, dir_fig):
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
    plt.xlim(0, Lseq + 10)
    plt.ylim(0, Lseq + 10)
    plt.title(f"<Neff>={average_neff}, <Neff>/L={average_neff / Lseq:.2f}")
    imgfile = os.path.join(dir_fig, f"afig1_top{top_2n}_{average_neff}_z{zcut}.png")
    plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def top_contact_map(dca_dataframe, pdbid, n_effective, xy_lim, pdb_dataframe, distance_cutoff, img_dir, theta):
    """
    Plots the top 10, 20, 30 DCA predictions and PDB monomer contact maps.
    :param img_dir:
    :param theta:
    :param dca_dataframe:
    :param pdbid:
    :param n_effective:
    :param xy_lim:
    :param pdb_dataframe:
    :param distance_cutoff:
    :return:
    """
    ii = 0
    n = 10
    plt.figure(ii + 1, figsize=(7, 7))
    colors = plt.cm.Paired(range(n + 1))
    # f = "smog.contacts"
    # df_smog = pd.read_csv(f, names=["chain_1", "i", "chain_2", "j"], delim_whitespace=True)
    plt.scatter('i', 'j', data=pdb_dataframe[pdb_dataframe["d"] <= distance_cutoff],
                color='lightgrey', s=85, label=f"{pdbid}")
    plt.scatter('j', 'i', data=pdb_dataframe[pdb_dataframe["d"] <= distance_cutoff],
                color='lightgrey', s=85, label=f"{pdbid}")
    # plt.scatter("j", "i", data=df_smog, s=85, label="smog")
    for jj in range(n, 0, -1):
        cn = int(jj * 10)
        plt.scatter('i', 'j', data=dca_dataframe[:cn], label=f"top{cn}", color=colors[jj],
                    edgecolors='black')
    plt.scatter('j', 'i', data=dca_dataframe[:300], label=f"top{300}", color='black',
                edgecolors='black')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=5, fancybox=True)
    plt.xlabel("resi")
    plt.ylabel("resj")
    plt.xlim(0, xy_lim + 2)
    plt.ylim(0, xy_lim + 2)
    if theta > 0:
        plt.title(f"theta={theta:.2f}")
        img_out = os.path.join(img_dir, f"top30_t{theta:.2f}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    else:
        img_out = os.path.join(img_dir, f"top30_neff{n_effective}.png")
        plt.title(f"Neff={n_effective}, Neff/L={n_effective / xy_lim:.2f}")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_score_distribution(dca_score, n_seqs, dca_dataframe, n_effective, n_res, img_dir, theta):
    plt.figure(n_seqs, figsize=(5,5))
    heights, bins, patches = plt.hist(dca_dataframe[dca_score], bins=50)
    plt.semilogy()
    if theta > 0:
        plt.title(f"theta={theta:.2f}")
    else:
        plt.title(f"Neff={n_effective}, Neff/L={n_effective / n_res:.2f}")
    plt.xlim(-10, 21)
    plt.xticks(np.arange(-10, 24, 3))
    mean = dca_dataframe[dca_score].mean()
    std = dca_dataframe[dca_score].std()
    plt.vlines(mean, 0, 1000, color="black", linestyles="dashed", label=f"mean={mean:.3f}, std={std:.3f}")
    plt.xlabel("Z-score")
    plt.ylabel("Counts")
    plt.legend(loc="upper right")
    if theta:
        img_out = os.path.join(img_dir, f"{dca_score}_theta{theta:.2f}.png")
        plt.savefig(img_out, format="png", dpi=200)
    else:
        img_out = os.path.join(img_dir, f"{dca_score}_neff{n_effective}.png")
        plt.savefig(img_out, format="png", dpi=200)
    plt.close()


def plot_dist_distribution(dca_dataframe, n_effective, n_res, img_dir):
    plt.figure(8927384929, figsize=(5,5))
    heights, bins, patches = plt.hist(dca_dataframe["d"], bins=50, color="red", edgecolor="black")
    plt.ylim(0, 20)
    plt.title(f"Neff={n_effective}, Neff/L={n_effective / n_res:.2f}")
    plt.xlabel("Intra-chain distance $\AA$")
    plt.ylabel("Counts")
    plt.legend(loc="best")
    img_out = os.path.join(img_dir, f"distance_neff{n_effective}.png")
    plt.savefig(img_out, format="png", dpi=200)
    plt.close()


def zscore_contact_map(dca_dataframe, pdbid, n_effective, xy_lim, pdb_dataframe, distance_cutoff,
                       img_dir, zcutoff, theta):
    """
    Plots the top 10, 20, 30 DCA predictions and PDB monomer contact maps.
    :param zcutoff:
    :param img_dir:
    :param theta:
    :param dca_dataframe:
    :param pdbid:
    :param n_effective:
    :param xy_lim:
    :param pdb_dataframe:
    :param distance_cutoff:
    :return:
    """
    i = 10000000
    plt.figure(i + 1, figsize=(7, 7))
    plt.scatter('i', 'j', data=pdb_dataframe[pdb_dataframe["d"] <= distance_cutoff],
                color='lightgrey', s=60, label=f"{pdbid}")
    plt.scatter('i', 'j', data=dca_dataframe, label=f"z={zcutoff}", color="blue",
                edgecolors='black')
    plt.legend(loc='best')
    plt.xlabel("resi")
    plt.ylabel("resj")
    plt.xlim(1, xy_lim + 2)
    plt.ylim(1, xy_lim + 2)
    if theta > 0:
        plt.title(f"theta={theta:.2f}")
        dir_cutoff = os.path.join(img_dir, f"z{zcutoff}")
        if not os.path.exists(dir_cutoff):
            os.makedirs(dir_cutoff)
        img_out = os.path.join(dir_cutoff, f"zscore_{zcutoff:.2f}_t{theta:.2f}.png")
        plt.savefig(img_out, format="png", dpi=200,
                    bbox_inches='tight')
    else:
        plt.title(f"Neff={n_effective}, Neff/L={n_effective / xy_lim:.2f}")
        img_out = os.path.join(img_dir, f"z{zcutoff}_neff{n_effective}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def plot_norms(a, pfamid, threshold_list, array_neff, dca_score, seq_len, zfilter, dir_avg_results):
    avg_neff = array_neff.mean(axis=0)
    neff_l_str = [str(f"{kk / seq_len:.3f}") for idk, kk in enumerate(avg_neff)]
    thresh_val_str = [str(f"{kk:.2f}") for idk, kk in enumerate(threshold_list)]
    # x_names = neff_l_str
    # hue_names = thresh_val_str
    # dim1, dim2, dim3 = np.meshgrid(x_names, np.arange(a.shape[1]), hue_names, indexing='ij')
    # sns.boxplot(x=dim1.ravel(), y=a.ravel(), hue=dim3.ravel())
    # plt.show()
    ma = np.mean(a, axis=1)
    sa = np.std(a, axis=1)
    max_ = np.max(a, axis=1)
    min_ = np.min(a, axis=1)
    plt.figure(0, figsize=(18, 10))
    for jj in range(len(ma)):
        # plt.errorbar(range(len(ma[jj])), ma[jj], yerr=max_[jj], capsize=6)
        # plt.errorbar(range(len(ma[jj])), ma[jj], yerr=min_[jj], capsize=6)
        plt.errorbar(range(len(ma[jj])), ma[jj], yerr=sa[jj], capsize=6)
        plt.scatter(range(len(ma[jj])), ma[jj], edgecolors="black", label=thresh_val_str[jj])

    plt.semilogy()
    plt.ylabel("Number of pairs")
    plt.xlabel("Neff/L")
    plt.legend(loc='best')
    plt.title(pfamid)
    plt.xticks(range(len(avg_neff)), neff_l_str)
    plt.grid(which="both", alpha=0.3, color="black")
    if zfilter:
        imgfile = os.path.join(dir_avg_results, f"normalization_zscore_ppv_100reps_{dca_score}.png")
        plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')


def average_ppv(pfam_id, metric_array, threshold_list, model_length_array,
                dca_score, seq_len, zfilter, dir_avg_results):
    # PLOT PPV vs. NEFF FOR DCA PREDICTIONS >= ZSCORE=5.6 and 3.5
    avg_neff = model_length_array.mean(axis=0)
    x = len(threshold_list)
    avg_ppv = np.mean(metric_array, axis=1)
    std_ppv = np.std(metric_array, axis=1)
    fig, ax = plt.subplots(figsize=(15, 7))
    # cmap = plt.get_cmap('viridis', 15)
    if zfilter:
        thresh_symbol = "z"
    else:
        thresh_symbol = "top_n"
    for ii in range(x):
        plt.errorbar(range(len(avg_ppv[ii])), avg_ppv[ii], yerr=std_ppv[ii], capsize=6)
        plt.scatter(range(len(avg_ppv[ii])), avg_ppv[ii], edgecolors="black",
                    label=f"{thresh_symbol}={threshold_list[ii]}")
    plt.ylabel("<PPV>")
    plt.xlabel("Neff/L")
    plt.legend(loc='best')
    plt.title(pfam_id)
    neff_l = [str(f"{kk / seq_len:.3f}") for idx, kk in enumerate(avg_neff)]
    plt.xticks(range(len(avg_neff)), neff_l)
    plt.grid(which="both", alpha=0.3, color="black")
    plt.ylim(0, 1)
    if zfilter:
        imgfile = os.path.join(dir_avg_results, f"zscore_ppv_100reps_{dca_score}.png")
        plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')
    else:
        imgfile = os.path.join(dir_avg_results, f"top{threshold_list[-1]}_ppv_100reps_{dca_score}.png")
        plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')
    plt.show()


def ranked_average_ppv(array_ppv, array_neff, dca_score, seq_len, dir_avg_results):
    avg_ppv = np.mean(array_ppv, axis=1)
    std_ppv = np.std(array_ppv, axis=1)
    # n_thresholds, n_models, n_pairs = avg_ppv.shape
    n_thresholds, n_models = avg_ppv.shape
    neff_l = [str(f"{kk / seq_len:.2f}") for idx, kk in enumerate(array_neff[0])]
    for i in range(n_models):
        neff_l = array_neff[0]
        plt.errorbar(range(len(avg_ppv[i])), avg_ppv[i], yerr=std_ppv[i],
                     capsize=6, label=f"$neff/L=${neff_l[i]}")
        plt.scatter(range(len(avg_ppv[i])), avg_ppv[i])
        # if float(neff_l[i]) < 1:
        #     plt.errorbar(range(len(avg_ppv[0][i])), avg_ppv[0][i], yerr=std_ppv[0][i],
        #                  capsize=6, linestyle="dashed", label=f"neff/L={neff_l[i]}")
        #     plt.scatter(range(len(avg_ppv[0][i])), avg_ppv[0][i])
        # else:
        #     plt.errorbar(range(len(avg_ppv[0][i])), avg_ppv[0][i], yerr=std_ppv[0][i],
        #                  capsize=6, label=f"$neff/L=${neff_l[i]}")
        #     plt.scatter(range(len(avg_ppv[0][i])), avg_ppv[0][i])
    plt.ylabel("<PPV>")
    plt.xlabel("ranked index")
    plt.ylim(0, 1)
    plt.legend(bbox_to_anchor=(1, 0.5), loc='center left')
    plt.grid(which="both", alpha=0.3, color="black")
    imgfile = os.path.join(dir_avg_results, f"ranked_ppv_100reps_{dca_score}.png")
    plt.savefig(imgfile, format="png", dpi=200, bbox_inches='tight')
    plt.show()


def plot_ppv(zfilter, ppv_list, neff_value, pfamid, sequence_len, threshold_val, img_dir):
    final_img_dir = img_dir + "\\ppv\\"
    if not os.path.exists(final_img_dir):
        os.makedirs(final_img_dir)
    plt.figure(figsize=(5, 5))
    plt.plot(ppv_list)
    plt.scatter(range(len(ppv_list)), ppv_list, edgecolors="black", label=f"{threshold_val:.1f}")
    plt.ylim(0, 1.1)
    plt.xlabel("ranked index")
    plt.ylabel("PPV")
    plt.title(f"{pfamid}, Neff:{neff_value}, Neff/L={neff_value / sequence_len:.2f}")
    plt.grid(which="both", alpha=0.3)
    plt.legend(loc="best")

    if zfilter:
        img_out = os.path.join(final_img_dir, f"ppv_neff{neff_value}_z{threshold_val:.1f}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    else:
        img_out = os.path.join(final_img_dir, f"ppv_neff{neff_value}_top{threshold_val:.1f}.png")
        plt.savefig(img_out, format="png", dpi=200, bbox_inches='tight')
    plt.close()


def pipeline_replicates(pfam_id, pdb_id, DIR_RESULTS, DIR_REPLICATES, DIR_PDB, dir_system, dca_score, ncols,
                        model_lengths,
                        neffective_array, thresholds_list, pdbid, seqlen, pfamid, pc, dca_mode,
                        zfilter=True, npairs=0, plots=False, load=False, passthrough=False):
    """
    For every replicate and model, process raw DCA results, calculate Zscores, and PPV
    :param dir_system:
    :param pfam_id:
    :param pdb_id:
    :param DIR_PDB:
    :param pc:
    :param dca_mode:
    :param DIR_RESULTS:
    :param DIR_REPLICATES:
    :param seqlen:
    :param pfamid:
    :param ncols:
    :param pdbid:
    :param npairs: to use in plotting ranked pairs
    :param load: Boolean; True then loads processed dca file.
    :param plots: Bool; Plots score distribution, top pairs, and zscore filtered pairs depending on zfilter.
    :param zfilter: Bool; If False, activates top pair program. If true, zscore filter program activates.
    :param passthrough: Bool; If True, skips this program and loads ppv array
    :param dca_score: Str; Score type.
    :param model_lengths: List; List of model lengths.
    :param neffective_array: Numpy Array; Multidimensional array of neff sequence lengths
    :param thresholds_list: List; List of Zscore cutoffs
    :return:
    """
    distance_cutoff = 8
    avg_neff = neffective_array.mean(axis=0)
    # CONTROL PPV OUTFILE NAME DEPENDING ON ZFILTER VALUE
    if zfilter:
        output_ppv = os.path.join(f"{DIR_RESULTS}\\average_ppv", f"{pdbid}_ppv_zscore_{dca_score}_100reps.npy")
        output_norm = os.path.join(f"{DIR_RESULTS}\\average_ppv", f"ppv_norm_z_{dca_score}_reps100.npy")
    else:
        if npairs > 0:
            output_ppv = os.path.join(f"{DIR_RESULTS}\\average_ppv",
                                      f"{pdbid}_ppv_top{npairs}_{dca_score}_100reps.npy")
        else:
            output_ppv = os.path.join(f"{DIR_RESULTS}\\average_ppv",
                                      f"{pdbid}_ppv_top{thresholds_list[-1]}_{dca_score}_100reps.npy")
            output_norm = os.path.join(f"{DIR_RESULTS}\\average_ppv",
                                       f"ppv_norm_top{thresholds_list[-1]}_{dca_score}_reps100.npy")

    # IF ALREADY CALCULATED LOAD PPV FILE
    if passthrough:
        return np.load(output_ppv), np.load(output_norm)

    # ELSE CALCULATE PPV
    else:
        n_replicates, n_sys = neffective_array.shape
        df_pdb = pipeline_pdb(pdbid, DIR_PDB)

        if npairs > 0:
            pos_pred_list = np.zeros((len(thresholds_list), n_replicates, n_sys, npairs))
            num_pairs_left = np.zeros_like(pos_pred_list)
        else:
            pos_pred_list = np.zeros((len(thresholds_list), n_replicates, n_sys))
            num_pairs_left = np.zeros_like(pos_pred_list)
        for rep_id in range(n_replicates):
            # Make directories for results and plots
            dir_dca_results = f"{DIR_REPLICATES}\\sub{rep_id}\\{dca_mode}\\pc{pc}\\"
            dir_contact_map = f"{dir_dca_results}images\\"
            if not os.path.exists(dir_contact_map):
                os.makedirs(dir_contact_map)

            for model_id in range(len(model_lengths)):
                n_degraded_seqs = avg_neff[model_id]
                raw_dca_output = os.path.join(dir_dca_results, f"DI_{pfamid}_n{n_degraded_seqs}.txt")

                neff = int(neffective_array[rep_id][model_id])
                # OUTPUT ZSCORE
                df_dca, map_idx = pipeline_process_results(pfam_id, pdb_id, raw_dca_output, dir_system, pc, df_pdb,
                                                           dca_score, neff, theta_flag=False, load_processed=load)
                draw_publish_dca(pfamid, df_dca, "A", dir_contact_map)
                # matrix = into_matrix(df_dca)
                assert df_dca is not None
                print(f"rep: {rep_id} model: {n_degraded_seqs} neff: {neff}  L: {ncols} map:{map_idx}")
                if plots:
                    # plot_score_distribution(dca_score, n_degraded_seqs, df_dca, neff, ncols, dir_contact_map)
                    if zfilter:
                        print("zfilter")
                        plot_dist_distribution(df_dca[:750], n_degraded_seqs, ncols, dir_contact_map)
                        plot_score_distribution("zscore", n_degraded_seqs, df_dca, neff, ncols,
                                                dir_contact_map, theta=False)
                        top_contact_map(df_dca, pdbid, neff, ncols + map_idx, df_pdb, distance_cutoff,
                                        dir_contact_map, theta=False)

                for k in range(len(thresholds_list)):
                    threshold_value = thresholds_list[k]
                    if zfilter:
                        df_filtered = df_dca[df_dca["zscore"] >= threshold_value]
                        if plots:
                            print(f"{threshold_value}")
                            # zscore_contact_map(df_filtered, pdbid, neff, ncols + map_idx, df_pdb, distance_cutoff,
                            #                    dir_contact_map, threshold_value, False)
                            # figure_1(df_dca, df_filtered, df_pdb, pdbid, distance_cutoff, n_degraded_seqs,
                            #          threshold_value, seqlen, map_idx, dir_contact_map)
                    else:
                        df_filtered = df_dca[:threshold_value]

                    if len(df_filtered) > 0:
                        # ppv, tpfp = calculate_ppv(df_filtered, distance_cutoff)
                        ppv=np.zeros(10)
                        tpfp=np.zeros(10)
                        if plots:
                            # plot_ppv(zfilter, ppv, neff, pfamid, seqlen, threshold_value, dir_contact_map)
                            print("PPV")
                        if npairs > 0:  # to use in plotting ranked pairs
                            pos_pred_list[k][rep_id][model_id] = ppv
                            num_pairs_left[k][rep_id][model_id] = tpfp
                        else:
                            pos_pred_list[k][rep_id][model_id] = ppv[-1]
                            num_pairs_left[k][rep_id][model_id] = tpfp[-1]
                    else:
                        pos_pred_list[k][rep_id][model_id] = 0
        # np.save(output_ppv, pos_pred_list)
        # np.save(output_norm, num_pairs_left)
        return pos_pred_list, num_pairs_left


def run_theta(pfam_id, DIR_RESULTS, msa_input, nthetas, pseudocount, passthrough=False):
    # THETA CHANGES
    dir_output_theta = os.path.join(f"{DIR_RESULTS}", "theta")
    neff_file = os.path.join(dir_output_theta, "neff_array_theta.npy")
    print(f"NEFF_FILE: {neff_file}")
    if passthrough:
        print(f"DIR: {dir_output_theta}")
        neff = np.load(neff_file)
    else:
        neff = np.zeros(nthetas)
        print("START MATLAB")
        for theta_idx in range(nthetas):
            theta = (theta_idx + 1) * 0.1
            if not os.path.exists(dir_output_theta):
                os.makedirs(dir_output_theta)
            dca_output = os.path.join(dir_output_theta, f"DI_{pfam_id}_t{theta:.2f}.txt")
            matrix_output = os.path.join(dir_output_theta, f"matrix_{pfam_id}_t{theta:.2f}.mat")
            if not os.path.exists(dca_output):
                print(f"theta:{theta}\tMSA_IN: {msa_input}\nDCA_OUT: {dca_output}\nMAT_OUT: {matrix_output}")
                eng = matlab.engine.start_matlab()
                eng.addpath(DIR_DCA_CODE, nargout=0)
                eng.addpath(dir_output_theta, nargout=0)
                neff[theta_idx] = eng.dca_pairdist(msa_input, dca_output, matrix_output, pseudocount,
                                                   theta, nargout=1)
            else:
                print(f"ALREADY RAN -> THETA:{theta}\tMSA_IN: {msa_input}\n"
                      f"DCA_OUT: {dca_output}\nMAT_OUT: {matrix_output}")
        print(f"DIR: {dir_output_theta}")
        np.save(neff_file, neff)
        return neff


def pipeline_reweighing(DIR_RESULTS, DIR_PDB, pfam_id, pdb_id, pc, dca_score, ncols, model_lengths, thresholds_list,
                        zfilter=True, npairs=0, plots=False, load=False, passthrough=False):
    distance_cutoff = 8
    DIR_THETA = os.path.join(DIR_RESULTS, "theta")
    DIR_CONTACT_MAP = os.path.join(DIR_THETA, "images")
    if not os.path.exists(DIR_CONTACT_MAP):
        os.makedirs(DIR_CONTACT_MAP)
    if zfilter:
        output_ppv = os.path.join(DIR_THETA, f"{pdb_id}_ppv_zscore_{dca_score}_theta.npy")
    else:
        if npairs > 0:
            output_ppv = os.path.join(DIR_THETA, f"{pdb_id}_ppv_top{npairs}_{dca_score}_theta.npy")
        else:
            output_ppv = os.path.join(DIR_THETA, f"{pdb_id}_ppv_top30_{dca_score}_theta.npy")

    if passthrough:
        return np.load(output_ppv)
    else:
        n_thetas = 10
        df_pdb = pipeline_pdb(pdb_id, DIR_PDB)

        if npairs > 0:
            pos_pred_list = np.zeros((len(thresholds_list), n_thetas, npairs))
        else:
            pos_pred_list = np.zeros(n_thetas)
            # pos_pred_list = np.zeros((len(thresholds_list), n_thetas))
        # Loop through theta values
        ppv = []
        for k in range(n_thetas):
            theta = (k + 1) * 0.1
            raw_dca_output = os.path.join(DIR_THETA, f"DI_{pfam_id}_t{theta:.2f}.txt")
            df_dca, map_idx = pipeline_process_results(pfam_id, pdb_id, raw_dca_output, DIR_SYS, pc, df_pdb, dca_score,
                                                       theta, theta_flag=True, load_processed=load)
            assert df_dca is not None
            print(f"theta: {theta}  L: {ncols}")
            if plots:
                # plot_score_distribution(dca_score, n_degraded_seqs, df_dca, neff, ncols, DIR_CONTACT_MAP)
                if zfilter:
                    plot_score_distribution("zscore", k, df_dca, 0, ncols, DIR_CONTACT_MAP, theta)
                    # top_contact_map(df_dca, pdb_id, 0, ncols + map_idx, df_pdb, distance_cutoff,
                    #                 DIR_CONTACT_MAP, theta)
            # Loop through every threshold and calculate PPV
            for ii in range(len(thresholds_list)):
                threshold_value = thresholds_list[ii]
                if zfilter:
                    df_filtered = df_dca[df_dca["zscore"] >= threshold_value]
                    if plots:
                        print(f"{threshold_value}")
                        # zscore_contact_map(df_filtered, pdb_id, 0, ncols + map_idx, df_pdb, distance_cutoff,
                        #                    DIR_CONTACT_MAP,
                        #                    threshold_value, theta)
                else:
                    df_filtered = df_dca[:threshold_value]
            # Calculate PPV and then save as npy file
            ppv.append([calculate_ppv(df_dca, distance_cutoff), df_dca["zscore"].to_numpy()])
            pos_pred = np.array(ppv)
        np.save(output_ppv, pos_pred)

    return pos_pred


def bootstrap(system_entry, dca_mode="mf", zbool=True, passthrough=False):
    pc = 0.2
    reweighting = 0.2
    nrep = 100
    score = "diapc"

    pfam_id = system_entry.pfam_id
    pdb_id = system_entry.pdb_id

    DIR_SYS, DIR_RESULTS, DIR_AVG_RESULTS, DIR_PDB, DIR_DIST_MAT, DIR_REPLICATES, nseq_list = make_dir_struct(pfam_id)
    file_raw_msa = os.path.join(DIR_SYS, f"{pfam_id}_full.txt")

    assert os.path.exists(file_raw_msa)

    file_filtered_msa = os.path.join(DIR_SYS, f"{pfam_id}_full_filtered_25.fasta")
    # 1. Gap filter
    if not os.path.exists(file_filtered_msa):
        file_filtered_msa, nseq, percent_gaps = gap_filter(file_raw_msa, 0.25)

    # Read in filtered MSA
    # 2. Check filesize
    seq_l = check_length(file_filtered_msa)

    # Read list of lengths generated by gen_replicates
    if os.path.exists(nseq_list):
        with open(nseq_list, "r") as fp:
            list_of_len = [int(line.rstrip()) for line in fp]
            assert list_of_len[-1] > 0
    else:
        # 3. Replicate generation
        list_of_len = generate_replicates(file_filtered_msa, nrep)

    # 4. RUN DCA FOR EACH REPLICATE ENSEMBLE
    if zbool:
        threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]
    else:
        threshold_values = np.arange(10, 110, 10)

    program = "replicates"
    if program == "replicates":
        neff_array = run_replicates(pfam_id, list_of_len, nrep, DIR_REPLICATES, pc, reweighting, passthrough=True)
        # neff_array = run_neff_calculation(list_of_len, nrep, DIR_REPLICATES, reweighting, passthrough=True)
        ppv_array, norm_array = pipeline_replicates(pfam_id, pdb_id, DIR_RESULTS, DIR_REPLICATES, DIR_PDB, DIR_SYS,
                                                    score, seq_l, list_of_len, neff_array, threshold_values, pdb_id,
                                                    seq_l, pfam_id, pc, dca_mode, zfilter=zbool, npairs=0, plots=True,
                                                    load=True, passthrough=False)
        if zbool:
            plot_norms(norm_array, pfam_id, threshold_values, neff_array, score, seq_l, zfilter=zbool,
                       dir_avg_results=DIR_AVG_RESULTS)
        average_ppv(pfam_id, ppv_array, threshold_values, neff_array, score, seq_l, zfilter=zbool,
                    dir_avg_results=DIR_AVG_RESULTS)
        # ranked_average_ppv(ppv_array, neff_array, score, seq_l)
    # if program == "theta":
    #                           run_theta(file_filtered_msa, 10, pc, passthrough=True)
    #                           ppv_array = pipeline_reweighing(score, seq_l, list_of_len, threshold_values, pdb_id, zfilter=zbool,
    #                           npairs=0, plots=False, load=True, passthrough=False)
    #                           n_thetas, nn = ppv_array.shape
    #                           for i in range(n_thetas):
    #                           theta_value = (i + 1) * 0.1
    #                           plt.plot(ppv_array[i][1], ppv_array[i][0])
    #                           plt.scatter(ppv_array[i][1], ppv_array[i][0], label=f"{theta_value:.2f}")
    #                           plt.legend(loc="best")
    #                           plt.ylim(0, 1)
    #                           plt.xlabel("zscore")
    #                           plt.ylabel("PPV")
    #                           plt.show()

    # plt.errorbar(range(10), avg_ppv[i], yerr=std_ppv[i], capsize=6)
    # plt.scatter(range(len(avg_ppv[i])), avg_ppv[i], edgecolors="black", label=f"top{z_cutoff_list[i]}")
    # plt.savefig(f"{DIR_RESULTS}theta\\zscore_ppv_theta.png", format="png", dpi=200, bbox_inches='tight')
    # D = defaultdict(list)
    # for i in range(10):
    #     for s in zl:
    #         D[f"{(i+1)*0.1}"].append(s["di"])


# INPUTS
def run_boot():
    sys_file = os.path.join(DIR_ROOT, "systems.csv")
    df_systems = pd.read_csv(sys_file, header=0)
    dca = "mf"
    # df_systems = df_systems[5:]
    for shp, sys_entry in df_systems.iterrows():
        print(f"{sys_entry.pfam_id}, {sys_entry.pdb_id}")
        if sys_entry.pfam_id == "PF00022":
            bootstrap(sys_entry, dca_mode=dca, zbool=True, passthrough=False)
        else:
            bootstrap(sys_entry, dca_mode=dca, zbool=True, passthrough=False)


run_boot()
