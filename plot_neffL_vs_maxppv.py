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
for k in range(len(threshold_values)):
    z = threshold_values[k]
    ppv = []
    neff_l = []
    avg_neff_l = []
    avg_ppv = []
    plt.figure(k + 1)
    for ix, i in enumerate(sysid[:num_sys]):
        ppv.append(np.load(f"results/{i}_ppv_zscore_diapc_100reps.npy"))
        neff_l.append(np.load(f"results/{i}_neff_array.npy"))
        avg_neff_l.append(np.mean(neff_l[ix][:, :] / l[ix], axis=0))
        avg_ppv.append(np.mean(ppv[ix][k, :, :], axis=0))
        # norm = 1 / np.max(avg_ppv[ix])
        # norm = 1
        # plt.plot(avg_neff_l[ix], avg_ppv[ix]*norm, label=f"{i}")
        # plt.scatter(avg_neff_l[ix], avg_ppv[ix]*norm)
        # plot max ppv vs neff/L
        ind = np.argmax(avg_ppv[ix])
        plt.scatter(avg_neff_l[ix][ind], np.max(avg_ppv[ix]), marker="s", label=f"{i}")

    plt.legend(loc="best")
    plt.title(f"Z={z}")
    plt.xlabel("Neff/L")
    plt.ylabel("max(PPV_avg)")
    plt.ylim(-0.1, 1.1)
    plt.semilogx()
    plt.xlim(-200, 200)
    plt.xticks([1e-2, 1e-1, 1e0, 1e1, 1e2])
    # plt.savefig(f"images/{num_sys}/z{z}_{num_sys}_normalized.png", format="png", dpi=150, bbox_inches='tight')
    plt.savefig(f"images/{num_sys}/maxppv_z{z}_{num_sys}.png")
    plt.close()
    # plt.show()
