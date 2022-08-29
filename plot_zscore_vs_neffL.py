import numpy as np
import matplotlib.pyplot as plt
from msa.tools.check_length import check_length

sysid = []
with open("sys_10.txt", "r") as fp:
    for line in fp.readlines():
        sysid.append(line.rstrip())
l= []
for jx, j in enumerate(sysid):
    l.append(check_length(f"aln/{j}.fa"))
threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]
num_sys = 10
zmax = np.zeros((len(threshold_values), num_sys))
for k in range(len(threshold_values)):
    z = threshold_values[k]
    ppv = []
    neff_l = []
    avg_neff_l = []
    avg_ppv = []
    # plt.figure(k + 1)
    for ix, i in enumerate(sysid[:num_sys]):
        ppv.append(np.load(f"results/{i}_ppv_zscore_diapc_100reps.npy"))
        neff_l.append(np.load(f"results/{i}_neff_array.npy"))
        avg_neff_l.append(np.mean(neff_l[ix][:, :] / l[ix], axis=0))
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
