import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from msa.tools.check_length import check_length
from analysis.plots import zscore_vs_distance
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

    zscore_vs_distance(df, system_id, neff, seq_l, path_results)
