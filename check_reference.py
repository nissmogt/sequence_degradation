import os
import pandas as pd
import analysis.zscore
import numpy as np
import matplotlib.pyplot as plt
from msa.tools.check_length import check_length

dca_dir = os.path.join("tests", "single")
sysid_list = []
with open("system_names.txt", "r") as fp:
    for line in fp.readlines():
        sysid_list.append(line.rstrip())

neff = []
scores = []
num_sys = len(sysid_list)
# num_sys = 76
systems = sysid_list[:num_sys]
seql_list = []
for jx, sysid in enumerate(systems):
    print(f"System ID: {sysid}")
    dir_sys = os.path.join(dca_dir, "systems", sysid)
    msa_file = os.path.join(dca_dir, "systems", sysid, f"{sysid}_filtered_25.fasta")
    # dir_results = os.path.join(os.path.dirname(msa_file), "results", "average_ppv")
    dir_results = os.path.join(os.path.dirname(msa_file), "results")
    # seql_list.append(check_length(msa_file))
    neff.append(np.load(os.path.join(dir_sys, "neff_array.npy")))
    seql_list.append(check_length(msa_file))
    if neff[jx] / seql_list[jx] > 1:
        df_filename = os.path.join(dir_sys, f"{sysid}_neff{neff[jx]}_pc0.2_all.txt")
        df_temp = pd.read_csv(df_filename, delim_whitespace=True)
        scores.append(df_temp["diapc"].to_numpy())

og_ref = np.load(os.path.join("assets", "monomer_DIAPC_scores_76.npy"))
ref_file = os.path.join("assets", "monomer_diapc.npy")
reference = analysis.zscore.make_reference(scores, "diapc", outf=ref_file)
# reference = np.load(ref_file)
plt.figure(0)
h1 = plt.hist(reference, density=True, bins=50,
              label=f"148 systems mean={np.mean(reference):.4f} std={np.std(reference):.4f}")
h2 = plt.hist(og_ref, density=True, bins=50, alpha=0.6, color="red",
              label=f"76 systems mean={np.mean(og_ref):.4f} std={np.std(og_ref):.4f}")
plt.semilogy()
plt.legend(loc="best")
plt.xlabel("diapc")
plt.ylabel("pdf")
plt.show()
