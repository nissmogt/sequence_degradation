import matplotlib.pyplot as plt
import numpy as np


def contact_map_single(couplings=None, monomer=None, n=10, x="text", distance_cutoff=6, symmetric=True):

    if len(couplings.columns) == 0 and len(monomer.columns) == 0:
        raise ValueError("Need to specify at least one of couplings or monomer")
    plt.style.use('./styles/contact_map_single.mplstyle')
    plt.figure(num=np.random.randint(10000))
    if len(monomer.columns) > 0:
        plt.scatter("i", "j", data=monomer[monomer["d"] <= distance_cutoff], label="monomer")
    plt.scatter("i", "j", data=couplings[:n], label=f"{x}")
    plt.legend(loc="best")
    plt.xlabel("resi")
    plt.ylabel("resj")
    plt.show()
