def plot_cm_mfdca(pdb_df_list, cutoff, length_a, length, atom, df_dca, msa_name=None,
                  other_dca=None, img_name=None, cbar=0):
    """
    :param df_dca:
    :param pdb_df_list:
    :param msa_name:
    :param cutoff:
    :param length_a:
    :param length:
    :param atom:
    :param other_dca:
    :return:
    """
    import matplotlib.pylab as plt
    import numpy as np
    from numpy.random import randint
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator
    plt.style.use('paper.mplstyle')

    print("\t\tPlotting pdb monomer and interface...")
    figNum = ord(msa_name[1]) * randint(100) + ord(msa_name[2]) * randint(101, 1000)  # figure number

    df_mon = pdb_df_list[0]
    df_inter = pdb_df_list[1]

    # Plotting
    fig, ax = plt.subplots(num=figNum)
    # monomer
    mon_cutoff = 12
    ax.scatter('i', 'j', data=df_mon[df_mon["d"] <= mon_cutoff],
               label='monomer {}$\AA$'.format(mon_cutoff),
               c='xkcd:grey', marker='o', s=50, alpha=0.7)
    # interface
    ax.scatter('i', 'j', data=df_inter[df_inter["d"] <= cutoff],
               label='dimer {}$\AA$'.format(cutoff),
               c='#CACACA', marker='s', s=50, alpha=0.7)

    if not other_dca.empty:
        ax.scatter('i', 'j', data=other_dca, label='top {}'.format(len(other_dca)),
                   # c='#322e2f', marker="s", s=150)
                   c='#0D1137', marker="s", s=150)

    if not df_dca.empty:
        import os
        label_ = (os.path.basename(img_name)).strip(msa_name)
        if cbar > 0:
            import matplotlib.colors
            norm = matplotlib.colors.Normalize(vmin=0.00284, vmax=1.1756)
            Reds = plt.get_cmap('Reds', 512)
            cmap = matplotlib.colors.ListedColormap(Reds(np.linspace(0.15, 1.3, 256)))
            sc = ax.scatter('i', 'j', data=df_dca, c=df_dca["di"].to_numpy(), cmap=cmap, vmin=0.00284, vmax=1.1756,
                            s=70, alpha=1, edgecolor='black', linewidth=0.2, label="DCA top {}".format(len(df_dca)))
            # color = 'orangered', s = 70)
            # ax.scatter('i', 'j', data=df_dca, label='top {}'.format(len(df_dca)),
            plt.colorbar(sc)
        else:
            ax.scatter('i', 'j', data=df_dca, label='Z-score {} top {}'.format(cbar, len(df_dca)), c='black', s=50)

    # Plot dimer separator line
    extend = 2
    # ax.hlines(length_a, 0, length, linestyles=(0, (5, 10)), alpha=0.7, colors='black')
    # ax.vlines(length_a, 0, length+extend, linestyles=(0, (5, 10)), alpha=0.7, colors='black')

    # plot design
    ax.legend(loc='best')
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
    if not df_dca.empty:
        imgname = "{}.svg".format(img_name)
    else:
        # Dont change this unless you want pdb distance matrix in a new directory
        img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\distance_matrix\\"
        imgname = "{}{}_len_{}_{}{}.tif".format(img_dir, msa_name, length, atom, cutoff)
    plt.savefig(imgname)
    # plt.show()
    plt.close()


