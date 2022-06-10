import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

DIR_ROOT = "D:\\bootstrap\\"
pfam_id = "PF00014"
L = 53
DIR_SYS = os.path.join(DIR_ROOT, pfam_id)
DIR_RESULTS = os.path.join(DIR_SYS, "results")
DIR_AVG_RESULTS = os.path.join(DIR_RESULTS, "average_ppv")
DIR_PDB = os.path.join(DIR_SYS, "PDB")
DIR_DIST_MAT = os.path.join(DIR_PDB, "distance_matrix")
DIR_REPLICATES = os.path.join(DIR_SYS, "replicates")
neff_file = os.path.join(DIR_REPLICATES, "neff_array.npy")
neff = np.load(neff_file)
neff_l = [str(f"{kk / L:.2f}") for idx, kk in enumerate(neff[0])]

ppv_1 = os.path.join(DIR_AVG_RESULTS, "5pti_ppv_zscore_diapc_100reps.npy")
ppv_2 = os.path.join(DIR_AVG_RESULTS, "5pti_ppv_top100_diapc_100reps.npy")
ppv_array1 = np.load(ppv_1)
ppv_array2 = np.load(ppv_2)
ap1 = ppv_array1.mean(axis=1)
ap2 = ppv_array2.mean(axis=1)

