def plot_cm(pdb_df_list, cutoff, length_a, length, atom, df_dca, msa_name=None):
    """
    :param df_dca:
    :param pdb_df_list:
    :param msa_name:
    :param cutoff:
    :param length_a:
    :param length:
    :param atom:
    :return:
    """
    import matplotlib.pylab as plt
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator
    # plt.style.use('paper.mplstyle')

    # Dont change this unless you want pdb distance matrix in a new directory
    img_dir = "images\\"
    imgname = "{}{}_len_{}_{}{}.tif".format(img_dir, msa_name, length, atom, cutoff)

    print("\t\tPlotting pdb monomer and interface...")
    figNum = 10

    # Separate monomer and interfacial contacts
    df_mon = pdb_df_list[0]
    df_inter = pdb_df_list[1]

    # Plotting
    fig, ax = plt.subplots(num=figNum)
    # monomer
    mon_cutoff = 8
    ax.scatter('i', 'j', data=df_mon[df_mon["d"] <= mon_cutoff],
               label='{} monomer, {}$\AA$'.format(msa_name[:4], mon_cutoff),
               c='xkcd:grey', marker='o', s=25, alpha=0.7)
    # interface
    ax.scatter('i', 'j', data=df_inter[df_inter["d"] <= cutoff],
               c='#CACACA', marker='s', s=25, alpha=0.7)
    if not df_dca.empty:
        ax.scatter('i', 'j', data=df_dca)

    extend = 2
    ax.hlines(length_a, 0, length, linestyles=(0, (5, 10)), alpha=0.7, colors='black')
    ax.vlines(length_a, 0, length+extend, linestyles=(0, (5, 10)), alpha=0.7, colors='black')

    # plot design
    # ax.legend(loc='best')
    ax.set_facecolor('#f5f0e1')
    # ax.set_facecolor('white')
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.major.formatter._useMathText = True
    plt.xlabel("monomer 1 residue index"), plt.ylabel("monomer 2 residue index")
    if atom == 'aa':
        plt.title("PDB: {} (aa)".format(msa_name[:4]))
    else:
        plt.title("PDB: {} (c$\alpha$)".format(msa_name))
    # plt.savefig(imgname, dpi=600)
    plt.show()
    # plt.close()


