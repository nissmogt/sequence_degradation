import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def run_difference(threshold, df1, df2):
    """ Establish an order parameter to quantify the differences between two contact maps
    """

    # Load contact maps
    # folder = "contact_maps\\"
    # # prefix = "atom_distance_matrix_"
    # prefix = "ca_distance_matrix_"
    # ext = ".pdb.txt"
    # pdb_id1 = "1L8C"
    # pdb_id2 = "1R8U"
    # f1 = "{}{}{}{}".format(folder, prefix, pdb_id1, ext)
    # f2 = "{}{}{}{}".format(folder, prefix, pdb_id2, ext)
    # threshold = 12

    df1 = df1[df1["chain_1"] == "A"]
    df1 = df1[df1["chain_1"] == df1["chain_2"]]
    df2 = df2[df2["chain_1"] == "B"]
    df2 = df2[df2["chain_1"] == df2["chain_2"]]
    df1 = df1[df1["d"] < threshold]
    df2 = df2[df2["d"] < threshold]
    df2.i = df2.i - 55
    df2.j = df2.j - 55

    # match pairs
    x = df1.merge(df2, on=["i", "j"], how="inner")
    d1 = x["d_x"].to_numpy()
    d2 = x["d_y"].to_numpy()
    pairs = x.iloc[:, :2].to_numpy()
    # pairs = pairs.transpose()

    # Quick check if residue names are the same
    for i, resi in enumerate(x.resnames_x):
        if resi != x.resnames_y[i]:
            print("False")

    # calculate difference in residue distances
    delta_d = d1 - d2
    p_array = np.column_stack((pairs, delta_d))
    df = pd.DataFrame(p_array, columns=["i", "j", "delta_d"])

    # plot heatmap
    fig, ax = plt.subplots(figsize=(5, 4))
    ax = plt.scatter("i", "j", data=df, c=df.delta_d, marker='s', s=8, cmap="jet")
    plt.colorbar()
    # plt.show()
    return fig
