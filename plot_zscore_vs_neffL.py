import numpy as np
import os
import matplotlib.pyplot as plt
from msa.tools.check_length import check_length

sysid = []
with open("system_names.txt", "r") as fp:
    for line in fp.readlines():
        sysid.append(line.rstrip())
seqlen_list= []
ppv_filenames = []
for jx, j in enumerate(sysid):
    seqlen_list.append(check_length(f"aln/{j}.fa"))
    ppv_filenames.append(os.path.join(f"results/{i}_ppv_zscore_diapc_100reps.npy"))
threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]
num_sys = len(sysid)
zmax = np.zeros((len(threshold_values), num_sys))
for k in range(len(threshold_values)):
    z = threshold_values[k]
    ppv = []
    neff_l = []
    avg_neff_l = []
    avg_ppv = []
    # plt.figure(k + 1)
    for ix, i in enumerate(sysid[:num_sys]):
        ppv.append(np.load())
        neff_l.append(np.load(f"results/{i}_neff_array.npy"))
        avg_neff_l.append(np.mean(neff_l[ix][:, :] / seqlen_list[ix], axis=0))
        avg_ppv.append(np.mean(ppv[ix][k, :, :], axis=0))
        ind = np.argmax(avg_ppv[ix])
        zmax[k, ix] = avg_neff_l[ix][ind]

plt.figure()
for kk in range(num_sys):
    plt.plot(threshold_values, zmax[:, kk])
    plt.scatter(threshold_values, zmax[:, kk], label=f"{sysid[kk]}")
plt.legend(loc="best")
plt.ylabel("Neff/L")
plt.xlabel("Z-score threshold")
plt.semilogy()
plt.yticks([1e-2, 1e-1, 1e0, 1e1, 1e2])
plt.show()
