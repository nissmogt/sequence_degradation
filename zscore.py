import os
import pandas as pd
import numpy as np
import time
import matplotlib.colors
import matplotlib.pyplot as plt

# from universal_reference_distribution import ref_dist

# 5 plots

rcutoff = 12


def plot_ppv_vs_zscore(df_ppv_list, labels, msaName, imgLocalDir, pairType):
    x1 = df_ppv_list[0]
    x2 = df_ppv_list[1]
    # plt.scatter(x1['zscore'], x1['ppv'], label=labels[0], color='orange')
    # plt.scatter(x2['zscore'], x2['ppv'], label=labels[1], marker='*')
    plt.scatter(x1['zscore'], x1['ppv'], label="seq/L >= 1 systems reference", color='orange')
    plt.scatter(x2['zscore'], x2['ppv'], label="individual sys reference", marker='*')
    plt.title(msaName)
    plt.ylabel('ppv')
    plt.xlabel('zscore')
    plt.legend(loc='best')
    imgName = "{}ppv_vs_zscore_{}_{}.png".format(imgLocalDir, pairType, msaName)
    plt.grid(axis='both', alpha=0.5)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close()


def plot_compare_zscore(dflist, msaName, labels, imgLocalDir, pairType):
    match_dfs = dflist[0].merge(dflist[1], how='inner', on=['i', 'j'], suffixes=('_score1', '_score2'))
    df_1 = match_dfs['zscore_score1']
    df_2 = match_dfs['zscore_score2']

    plt.scatter(df_1, df_2, alpha=0.7)
    maxzscore = int(max(df_1.max(), df_2.max()))
    minzscore = int(min(df_1.min(), df_2.min()))
    plt.plot(range(minzscore, maxzscore), range(minzscore, maxzscore), color='black')
    plt.grid(which='both', alpha=0.2)
    plt.xlabel("seq/L >= 1 systems reference")
    plt.ylabel("individual sys reference")
    # plt.xlabel(labels[0])
    # plt.ylabel(labels[1])
    plt.title("{} {}".format(msaName, pairType))
    imgName = "{}compare_zscores_{}_{}.png".format(imgLocalDir, pairType, msaName)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close()


def plot_npairs_vs_zscore(df, labels, msa, i, imgLocalDir, pairType):
    plt.figure(100000 * (i + 1), figsize=(10, 5), dpi=100)
    plt.scatter(range(len(df[0])), y='zscore', data=df[0], c='orange', label="seq/L >= 1 systems reference")
    plt.scatter(range(len(df[1])), y='zscore', data=df[1], marker='.', label="individual sys reference", alpha=0.7)
    # plt.scatter(range(len(df[0])), y='zscore', data=df[0], c='orange', label=labels[0])
    # plt.scatter(range(len(df[1])), y='zscore', data=df[1], marker='.', label=labels[1], alpha=0.7)
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("z-score")
    plt.title("{}".format(msa))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    # plt.show()
    imgName = "{}zscore_vs_npairs_{}_{}.png".format(imgLocalDir, pairType, msa)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close(100000 * (i + 1))


def plot_histogram_scores(dataDf, labelList, msaName, idx, imgLocalDir, pairType):
    _x = dataDf[0]['zscore']
    _y = dataDf[1]['zscore']
    plt.figure(20000 * (idx + 1), figsize=(10, 5), dpi=100)
    bins = np.histogram(np.hstack((_x, _y)), bins=50)[1]
    plt.hist(_x, bins=bins, edgecolor='black', label="seq/L >= 1 systems reference", color='orange')
    plt.hist(_y, bins=bins, edgecolor='black', label="individual sys reference", alpha=0.5)
    # plt.hist(_x, bins=bins, edgecolor='black', label=labelList[0], color='orange')
    # plt.hist(_y, bins=bins, edgecolor='black', label=labelList[1], alpha=0.5)
    plt.yscale('log')
    plt.legend(loc='best')
    plt.xlabel('zscore')
    plt.ylabel("counts")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgName = "{}hist_{}_{}.png".format(imgLocalDir, pairType, msaName)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close(20000 * (idx + 1))


def plot_distance_vs_zscore(df, labels, msa, i, imgLocalDir, pairType):
    plt.figure(10000 * (i + 1), figsize=(10, 5), dpi=100)
    # plt.scatter(x='d', y='zscore', data=df[0], c='orange', label=labels[0])
    # plt.scatter(x='d', y='zscore', data=df[1], marker='.', label=labels[1], alpha=0.7)
    plt.scatter(x='d', y='zscore', data=df[0], c='orange', label="seq/L >= 1 systems reference")
    plt.scatter(x='d', y='zscore', data=df[1], marker='.', label="individual sys reference", alpha=0.7)
    plt.xlabel("distance")
    plt.ylabel("z-score")
    plt.title("{}".format(msa))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    # plt.show()
    imgName = "{}distance_vs_zscore_{}_{}.png".format(imgLocalDir, pairType, msa)
    # imgName = "{}{}_noAPC_monomer_distance_vs_zscore_vanilla_vs_nonbonded_restraint12A.png".format(imgLocalDir, msa)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close(10000 * (i + 1))


def plot_std(std_list, labelList, msaName, idx, scoreType, imgLocalDir):
    plt.figure(100 * (idx + 1), figsize=(10, 5), dpi=100)
    plt.plot(range(len(std_list[0])), std_list[0], c='orange', label=labelList[0])
    plt.plot(range(len(std_list[1])), std_list[1], label=labelList[1])
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("standard deviation (zscore)")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgName = "{}{}_removed_dca_vs_std_{}_{}.png".format(imgLocalDir, msaName, pair_type, scoreType)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close(100 * (idx + 1))


def plot_mean(mean_list, labelList, msaName, idx, imgLocalDir):
    plt.figure(2000 * (idx + 1), figsize=(10, 5), dpi=100)
    plt.plot(range(len(mean_list[0])), mean_list[0], c='orange', label=labelList[0])
    plt.plot(range(len(mean_list[1])), mean_list[1], label=labelList[1])
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("mean (zscore)")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgName = "{}{}_removed_dca_vs_mean_{}.png".format(imgLocalDir, msaName, pair_type)
    plt.savefig(imgName, dpi=900, bbox_inches='tight')
    plt.close(2000 * (idx + 1))


def stats_removed_pairs(df):
    _std = np.zeros(len(df))
    _mean = np.zeros(len(df))
    print(len(df))
    for p in range(len(_std)):
        _std[p] = df['zscore'][p:].std()
        _mean[p] = df['zscore'][p:].mean()
    return _std, _mean


def zscore_calc(data, ref, which, score, stats=False):
    if which == 'monomer':
        _df = data[data['chain_1'] == data['chain_2']].reset_index(drop=True)
    if which == 'interface':
        _df = data[data['chain_1'] != data['chain_2']].reset_index(drop=True)
    if which == 'all':
        _df = data
    if len(ref) > 0:
        _df['zscore'] = (_df[score] - ref.mean()) / ref.std()
    else:
        # _df['zscore'] = (_df[score] - _df[score].mean()) / _df[score].std()
        _df['zscore'] = (_df[score] - data[score].mean()) / data[score].std()

    if stats:
        _s, _m = stats_removed_pairs(_df)
        return _df, _s, _m
    else:
        return _df


def test1():
    from universal_reference_distribution import ref_dist
    home_dir = "D:\\dca-interface_copy\\"
    results_directory = ["nonbonded_restraints_results\\12A\\", "sasa_nonbonded_restraints_results\\",
                         "sasa_restraints_results\\sasa_5\\", "mfdca_results\\", "vanilla_results\\"]

    dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
              "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A',
              '5M72_A_5M72_B', '1EFP_A_1EFP_B', '1EP3_A_1EP3_B', '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B',
              '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']
    rd1 = home_dir + results_directory[4]
    df = pd.read_csv(rd1 + dimers[1], delimiter="\t", header=0)
    z = zscore_calc(df, ref=ref_dist(), which="interface")


def test():
    home_dir = "D:\\dca-interface_copy\\"
    results_directory = ["nonbonded_restraints_results\\12A\\", "sasa_nonbonded_restraints_results\\",
                         "sasa_restraints_results\\sasa_5\\", "mfdca_results\\", "vanilla_results\\"]
    # results_directory = ["vanilla_results\\", "sasa_restraints_results\\sasa_5\\"]

    dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
              "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A',
              '5M72_A_5M72_B', '1EFP_A_1EFP_B', '1EP3_A_1EP3_B', '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B',
              '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']
    # dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
    #           "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
    #           '5MU7_B_5MU7_A', '5M72_A_5M72_B',  '1EFP_A_1EFP_B', '1EP3_A_1EP3_B',
    #           '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B', '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']
    # dimers = ['1EFP_A_1EFP_B', '1EP3_A_1EP3_B', '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B',
    #           '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']
    # dimers = ['3MML_A_3MML_B', '1EP3_A_1EP3_B', '5L8H_B_5L8H_A', '1ZUN_A_1ZUN_B']
    rd1 = home_dir + results_directory[4]
    rd2 = home_dir + results_directory[4]

    l = [rd1.split('\\')[0], rd2.split('\\')[0]]
    print(rd1, rd2)

    reweights = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # pt = ['monomer', 'interface']
    pt = ['interface']
    sys_idx = 1
    start_time = time.time()
    for sys_idx in range(len(dimers)):
        score = 'fn_apc'
        # for idx, k in enumerate(reweights):
        for pair_type in pt:
            nPairs = 10
            sys_name = dimers[sys_idx]
            img = "D:\\dca-interface_copy\\test\\"
            img_local_dir = "{}{}\\{}\\{}\\{}\\ZSCORE_THRESH\\".format(img, sys_name[:4], l[1], pair_type, score.upper())
            if not os.path.exists(img_local_dir):
                os.makedirs(img_local_dir)
            print(sys_name, pair_type)
            data = "{}{}_{}_{}_mapped_aa_dist_zscore_ref_all.txt".format(rd1, score.upper(), pair_type, sys_name)
            data2 = "{}{}_{}_{}_mapped_aa_dist_zscore.txt".format(rd2, score.upper(), pair_type, sys_name)
            df1 = pd.read_csv(data, delimiter='\t', header=0)
            df2 = pd.read_csv(data2, delimiter='\t', header=0)
            df1 = df1[df1["zscore"] >= 3][:10]
            df2 = df2[df2["zscore"] >= 3][:10]

            df_list = [df1[:], df2[:]]

            if pair_type == 'monomer':
                threshold = 8
            else:
                threshold = 12
            df_ppv = ppv_calc(df_list[0], threshold)
            df_ppv2 = ppv_calc(df_list[1], threshold)

            # x1 = np.array([df1['zscore'].to_numpy(), df_ppv])
            # x2 = np.array([df2['zscore'].to_numpy(), df_ppv2])
            x = [df_ppv, df_ppv2]
            plot_ppv_vs_zscore(x, l, sys_name, img_local_dir, pair_type)

            # plt.plot(df_ppv.loc[:, 'ppv'], color='orange', label=l[0], marker='o')
            # plt.plot(df_ppv2.loc[:, 'ppv'], label=l[1], marker='*')
            plt.plot(df_ppv.loc[:, 'ppv'], color='orange', label="seq/L >= 1 systems reference", marker='o')
            plt.plot(df_ppv2.loc[:, 'ppv'], label="individual system reference", marker='*')
            plt.title(sys_name)
            plt.xlabel("rank-order pairs")
            plt.ylabel("ppv")
            # plt.xscale('log')
            plt.legend(loc='best')
            plt.grid(axis='both', alpha=0.5)
            imgName = "{}ppv_vs_npairs_{}_{}.png".format(img_local_dir, pair_type, sys_name)
            plt.savefig(imgName, dpi=900, bbox_inches='tight')
            plt.close()

            plot_compare_zscore(df_list, sys_name, l, img_local_dir, pair_type)
            plot_distance_vs_zscore(df_list, l, sys_name, sys_idx, img_local_dir, pair_type)
            plot_npairs_vs_zscore(df_list, l, sys_name, sys_idx, img_local_dir, pair_type)
            plot_histogram_scores(df_list, l, sys_name, sys_idx, img_local_dir, pair_type)

            # if score_type == 'fn_apc':
            # header = "i\tj\tfn_apc\tfn\tui\tuj\td\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id\tzscore"
            # out1 = "{}FN_{}_{}_mapped_aa_dist_zscore.txt".format(rd1, pair_type, sys_name)
            # out2 = "{}FN_{}_{}_mapped_aa_dist_zscore.txt".format(rd2, pair_type, sys_name)
            # df1.to_csv(out1, sep='\t', index=False, header=header, float_format='%.5f')
            # df2.to_csv(out2, sep='\t', index=False, header=header, float_format='%.5f')

    print("----elapsed time {}".format(time.time() - start_time))