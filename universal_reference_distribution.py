def ref_dist(score='fn_apc'):
    """
    Reference distribution for calculating Z-score based on DCA score.
    :param rd: String - DCA results directory
    :param ptype: String - Type of reference: {'interface', 'all'}
    :param score: String - Score definition from DCA output

    """
    # score = 'fn_apc'
    import pandas as pd
    import numpy as np
    import os
    file = "dimers_stats.csv"
    sys_df = pd.read_csv(file, header=0)
    high_med_df = sys_df.loc[:, ("MSA", "Neff/L")][sys_df["Neff/L"] >= 1]["MSA"]
    # high_df = sys_df.loc[:, ("MSA", "Neff/L")][sys_df["Neff/L"] >= 10]["MSA"]
    # t1 = sys_df.loc[:, ("MSA", "Neff/L")][sys_df["Neff/L"] < 10]
    # med_df = t1.loc[:, ("MSA", "Neff/L")][t1["Neff/L"] >= 1]["MSA"]
    # low_df = sys_df.loc[:, ("MSA", "Neff/L")][sys_df["Neff/L"] < 1]["MSA"]

    high_med = high_med_df.to_list()
    # high = high_df.to_list()
    # med = med_df.to_list()
    # low = low_df.to_list()
    dimers = sys_df["MSA"].to_list()
    # using dimers ref MEAN:2.231226761403696e-05    STD: 0.024543407209537002
    scores = []
    # rd = "D:\\dca-interface_copy\\vanilla_results\\"
    rd = "D:\\dca-interface_copy\\mfdca_results\\"
    lsys = high_med
    for i in range(len(lsys)):
        sys_name = lsys[i]

        # data = "{}{}_{}_mapped_aa_dist.txt".format(rd, score.upper(), sys_name)
        data = "{}DI_{}_mapped_aa_dist.txt".format(rd, sys_name)
        if os.path.exists(data):
            df = pd.read_csv(data, delimiter='\t', header=0)
            x = df[df["chain_1"] == df["chain_2"]]
            scores.append(x.loc[:, score].to_numpy())
    reference = np.concatenate(scores)
    np.save(f"monomer_{score.upper()}_scores.npy", reference)
    return reference


ref_dist("di_apc")

# grain = 'aa'
# p1 = "structures\\1R8U.pdb"
# p2 = "structures\\1L8C.pdb"
# d1 = "FN_APC_TAZ1_CITED2_hmmer.fas.txt"
# d2 = "FN_APC_TAZ1_HIF1A_hmmer.fas.txt"
# f1 = "1R8U_all_dca.txt"
# f2 = "1L8C_all_dca.txt"

# f1 = "{}_{}_inter_mapped_{}_dist.txt".format(os.path.basename(p1).strip(".pdb"), d1.strip(".fas.txt"), grain)
# f2 = "{}_{}_inter_mapped_{}_dist.txt".format(os.path.basename(p2).strip(".pdb"), d2.strip(".fas.txt"), grain)

# PLOT HISTOGRAM
# import matplotlib.pyplot as plt
# scores_cat1 = np.concatenate(scores)
# d = [f1, f2]
# s = []
# for i in range(2):
#     df = pd.read_csv(d[i], delimiter='\t', header=0)
#     x = df
# s.append(x.loc[:, score].to_numpy())
# scores.append(x.loc[:, score].to_numpy())
# scores_cat2 = np.concatenate(scores)
# bins = np.linspace(min(scores_cat2), max(scores_cat2), 100)
# b = plt.hist(scores_cat1, bins=bins, label="Original Zscore Reference")
# b2 = plt.hist(scores_cat2, bins=bins, label="Zscore Reference w/ 1R8U and 1L8C systems", alpha=0.6)
# plt.semilogy()
# plt.xlabel("CN-score")
# plt.ylabel("Counts")
# plt.legend(loc="best")
# plt.show()

# scores_cat = np.concatenate(scores)
# return scores_cat
