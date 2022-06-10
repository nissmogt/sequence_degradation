import pandas as pd


def draw_publish_dca(msaName, df_dca, chains, resDir):
    """
    Draws pairs on VMD structure.
    Author: Kareem Mehrabiani
    Date: March 3rd, 2020

    Updated to version 1.0 - October 21, 2020

    Usage: draw_contacts.py -h for help
    """
    # pdb_file = "distance_matrix\\heavy_atom_distance_matrix_{}.txt".format(msaName)
    # df_pdb = pd.read_csv(pdb_file, delimiter='\t', usecols=(0, 1, 3, 4))

    # vmd_mappings = df_dca.merge(df_pdb, how='inner', on=['i', 'j'])

    nPairs = len(df_dca)
    import os
    outDir = os.path.join(resDir, "tcl_script")

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    tcl_file = "{}draw_{}_top{}.tcl".format(outDir, msaName, nPairs)
    target = open("{}".format(tcl_file), 'w', encoding='utf-8')
    target.write("color Display Background white\n")
    target.write("axes location Off\n")
    target.write("display projection Perspective\n")
    target.write("display ambientocclusion on\n")
    target.write("display depthcue on\n")
    target.write("display rendermode GLSL\n")
    target.write("mol material AOShiny\n")
    atomselect = 0
    for i in range(nPairs):
        resi, resj = df_dca["si"][i], df_dca["sj"][i]
        # atomi, atomj = eval(df_dca["atom_id"][i])[0], eval(df_dca["atom_id"][i])[1]

        # residue i
        target.write("mol representation NewCartoon\n")
        target.write("mol material AOShiny\n")
        target.write("mol color ColorID 30\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resi, chains[0]))
        target.write("mol addrep top\n")
        target.write("mol representation Surf\n")
        target.write("mol material AOChalky\n")
        target.write("mol color ColorID 30\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resi, chains[0]))
        target.write("mol addrep top\n")
        target.write("set sel%d%d%s [atomselect top \"resid %d and chain %s\"]\n" % (i, resi, chains[0],
                                                                                     resi, chains[0]))
        target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
        atomselect += 1
        # residue j
        target.write("mol representation NewCartoon\n")
        target.write("mol material AOShiny\n")
        target.write("mol color ColorID 19\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resj, chains[1]))
        target.write("mol addrep top\n")
        target.write("mol representation Surf\n")
        target.write("mol material AOChalky\n")
        target.write("mol color ColorID 15\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resj, chains[1]))
        target.write("mol addrep top\n")
        target.write("set sel%d%d%s [atomselect top \"resid %d and chain %s \"]\n" % (i, resj, chains[1],
                                                                                      resj, chains[1]))
        target.write("lassign [atomselect%d get {x y z}] pos2\n" % atomselect)
        atomselect += 1

        target.write("draw color cyan\n")
        target.write("draw line $pos1 $pos2 style dashed width 5\n")

    # Two chains main representation and color
    target.write("mol modselect 0 top \"chain %s\"\n" % chains[0])
    target.write("mol modmaterial 0 top AOEdgy\n")
    target.write("mol modstyle 0 top NewCartoon\n")
    target.write("mol modcolor 0 top ColorID 6\n")
    target.write("mol modselect 1 top \"chain %s\"\n" % chains[1])
    target.write("mol modmaterial 1 top AOEdgy\n")
    target.write("mol modstyle 1 top NewCartoon\n")
    target.write("mol modcolor 1 top ColorID 9\n")
    target.write("display resetview\n")
    target.close()


def draw_dca(msaName, df_dca, chains, resDir):
    """
    Draws pairs on VMD structure.
    Author: Kareem Mehrabiani
    Date: March 3rd, 2020

    Updated to version 1.0 - October 21, 2020

    Usage: draw_contacts.py -h for help
    """
    # pdb_file = "distance_matrix\\heavy_atom_distance_matrix_{}.txt".format(msaName)
    # df_pdb = pd.read_csv(pdb_file, delimiter='\t', usecols=(0, 1, 3, 4))

    # vmd_mappings = df_dca.merge(df_pdb, how='inner', on=['i', 'j'])

    nPairs = len(df_dca)
    import os
    outDir = "tcl_scripts\\vanilla_results\\"
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    tcl_file = "{}draw_{}_FN_APC_top{}.tcl".format(outDir, msaName, nPairs)
    target = open("{}".format(tcl_file), 'w', encoding='utf-8')
    target.write("color Display Background white\n")
    target.write("axes location Off\n")
    target.write("display projection Orthographic\n")
    target.write("display ambientocclusion on\n")
    target.write("display depthcue on\n")
    target.write("display rendermode GLSL\n")
    target.write("mol representation VDW\n")
    target.write("mol material AOChalky\n")
    atomselect = 0
    for i in range(nPairs):
        resi, resj = df_dca["si"][i], df_dca["sj"][i]
        atomi, atomj = eval(df_dca["atom_id"][i])[0], eval(df_dca["atom_id"][i])[1]

        # residue i
        target.write("mol color ColorID 31\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resi, chains[0]))
        target.write("mol addrep top\n")
        target.write("mol color Element\n")
        target.write("mol selection {resid %d and chain %s and name %s}\n" % (resi, chains[0], atomi))
        target.write("mol addrep top\n")
        target.write("set sel%d%d%s [atomselect top \"resid %d and chain %s and name %s\"]\n" % (i, resi, chains[0],
                                                                                                 resi, chains[0],
                                                                                                 atomi))
        target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
        atomselect += 1
        # residue j
        target.write("mol color ColorID 15\n")
        target.write("mol selection {resid %d and chain %s}\n" % (resj, chains[1]))
        target.write("mol addrep top\n")
        target.write("mol color Element\n")
        target.write("mol selection {resid %d and chain %s and name %s}\n" % (resj, chains[1], atomj))
        target.write("mol addrep top\n")
        target.write("set sel%d%d%s [atomselect top \"resid %d and chain %s and name %s\"]\n" % (i, resj, chains[1],
                                                                                                 resj, chains[1],
                                                                                                 atomj))
        target.write("lassign [atomselect%d get {x y z}] pos2\n" % atomselect)
        atomselect += 1

        target.write("draw color green\n")
        target.write("draw line $pos1 $pos2 style dashed width 4\n")

    # Two chains main representation and color
    target.write("mol modselect 0 top \"chain %s\"\n" % chains[0])
    target.write("mol modmaterial 0 top AOEdgy\n")
    target.write("mol modstyle 0 top NewCartoon\n")
    target.write("mol modcolor 0 top ColorID 6\n")
    target.write("mol modselect 1 top \"chain %s\"\n" % chains[1])
    target.write("mol modmaterial 1 top AOEdgy\n")
    target.write("mol modstyle 1 top NewCartoon\n")
    target.write("mol modcolor 1 top ColorID 9\n")
    target.write("display resetview\n")
    target.close()


def batch_draw(f_list_msa, resDir):
    for i in range(len(f_list_msa)):
        msa_name = f_list_msa[i]
        chainList = [msa_name.split("_")[1], msa_name.split("_")[3]]
        data = "{}FN_APC_interface_{}_mapped_aa_dist_zscore_ref_highmed_all.txt".format(resDir, msa_name)
        df_inter = pd.read_csv(data, delimiter='\t', header=0)
        df_inter = df_inter[6].reset_index(drop=True)
        draw_publish_dca(msa_name, df_inter, chainList, resDir)


# Test
# dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
#           "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A',
#           '5M72_A_5M72_B', '1EFP_A_1EFP_B', '1EP3_A_1EP3_B', '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B',
#           '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']
# rdir = "vanilla_results\\"
# file = "dimers_stats.csv"
# sys_df = pd.read_csv(file, header=0)
# dimers = sys_df["MSA"].to_list()
# batch_draw(dimers, rdir)
# data = "{}{}_neff542_pc0.2_all.txt".format(resDir, pfam_id)
# df_inter = pd.read_csv(data, delimiter='\t', header=0)

# import pandas as pd
# msaname = "3A0R_A_3A0R_B"
# ch = [msaname.split("_")[1], msaname.split("_")[3]]
# data = "{}FN_APC_interface_{}_mapped_aa_dist_zscore_ref_all.txt".format(rdir, msaname)
# df = pd.read_csv(data, delimiter='\t', header=0)
# df = df[df["zscore"] >= 4.2][:5]
# draw_dca(msaname, df, ch, rdir)
