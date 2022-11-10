import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from msa.tools.check_length import check_length
import os

GOOGLEDRIVE_PATH = os.path.join("/Volumes", "GoogleDrive", "My Drive", "sequence_degradation")
sysid_list = []
with open("system_names.txt", "r") as fp:
    for line in fp.readlines():
        sysid_list.append(line.rstrip())

for system_id in sysid_list[:10]:
# system_id = "1cc8A"
    print(system_id)
    path_results = os.path.join(GOOGLEDRIVE_PATH, "systems", system_id)
    path_old_results = os.path.join(GOOGLEDRIVE_PATH, "systems_old_ref", system_id)

    if not os.path.exists(path_results):
        print(f"File {path_results} doesn't exist.")

    neff = np.load(os.path.join(path_results, "neff_array.npy"))
    msa_file = os.path.join(path_results, f"{system_id}_filtered_25.fasta")
    seq_l = check_length(msa_file)
    dca_file = os.path.join(path_results, f"{system_id}_neff{neff}_pc0.2_all.txt")
    dca_file_compare = os.path.join(path_old_results, f"{system_id}_neff{neff}_pc0.2_all.txt")
    df = pd.read_csv(dca_file, delimiter="\t")
    df_compare = pd.read_csv(dca_file_compare, delimiter="\t")

    plt.scatter(df["zscore"], df["d"], color="indigo", label="148 psicov monomers", marker="s")
    # plt.scatter(df_compare["zscore"], df_compare["d"], color="yellow", edgecolors="black", alpha=0.7,
    #             label="76 evcoupling (monomers)")
    # max_z = max(max(df["zscore"]), max(df_compare["zscore"]))
    max_z = max(df["zscore"])
    plt.hlines(8, 0, max_z, colors="xkcd:green", linestyles="dashed")
    plt.hlines(10, 0, max_z, colors="xkcd:black", linestyles="dashed")
    plt.hlines(12, 0, max_z, colors="xkcd:red", linestyles="dashed")
    plt.legend(loc="best")
    plt.xlabel("zscore")
    plt.ylabel("distance (A)")
    plt.title(f"{system_id} Neff/L={neff/seq_l:.2f}")
    # img_path = os.path.join(f"{path_results}", "images", f"compare_z_distance_evfold_psicov_systems.png")
    img_path = os.path.join(f"{path_results}", "images", f"z_vs_distance.png")
    plt.savefig(img_path, format="png", dpi=150,
                bbox_inches='tight')
    # plt.show()
    plt.close()
