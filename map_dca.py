#!/home/kmm5/anaconda2/bin/python
##############################################################################################
#
# Author: Faruck Morcos
# Revised on 10/02/13 by Kareem Mehrabiani 
#     -- added argparse to make userx easier
# Usage: python map_dca.py [hmmscan file] [dca contact map]
# Email: kmm5@rice.edu for questions or improvements.        
#
##############################################################################################


# take the name of files
# parser = argparse.ArgumentParser(description="Maps Protein Domain to a PDB structure")
# parser.add_argument("hmmscan_file", help="Hmmscan output")
# parser.add_argument("dca_pairs", help="DCA pairs is a three-column file with i and j\
#         equal to residue i and j and the third column is the DCA score.")
# args = parser.parse_args()

# align = args.hmmscan_file
# ranked = args.dca_pairs


def mapit(dca_data, pdbid, outdir, align):
    import linecache
    import numpy as np
    import pandas as pd
    import os

    # align = "2zwh_scan.txt"
    if not os.path.exists(align):
        print(f"{align} doesnt exist!")

    # get position of alignment data
    w = 1

    while len(linecache.getline(align, w).split("== domain")) < 2:
        w += 1

    # get information from alignment file
    ##NOTE: sometimes hmmscan output format is displaced by 1 row. So if there are any errors
    #      displayed in the output of this script, change r to 2.
    map_ref = os.path.join(outdir, f"{pdbid}_map_reference.txt")
    if os.path.exists(map_ref):
        df = pd.read_csv(map_ref, header=0, delimiter="\t")
        domain_aacode = df["Domain"].to_numpy()
        domain_num = df["n"].to_numpy()
        protein_aacode = df["Protein"].to_numpy()
        protein_num = df["n.1"].to_numpy()

    else:
        disp = 1
        domain = linecache.getline(align, w + disp).split()
        prot = linecache.getline(align, w + disp + 2).split()

        print('\n\tDomain:\t' + domain[0].rjust(10) + '\n\tProtein ID:\t' + prot[0].rjust(10) + '\n')
        print('\tGenerating matched alignment...\n')

        d_init = int(domain[1])
        d_end = int(domain[3])
        p_init = int(prot[1])
        p_end = int(prot[3])
        # print("%s %s %s %s\n" % (d_init, d_end-d_init, p_init, p_end-p_init))
        # domain sequence string
        l1 = domain[2]
        # protein sequence string
        l2 = prot[2]

        x1 = len(l1)
        x2 = len(l2)

        output1 = open(map_ref, "w")
        output1.write('Domain\tn\tProtein\tn\n')

        # get the difference between initial positions
        # delta = max(x1,x2)
        delta = max(d_end - d_init, p_end - p_init)

        # domain code and respective number
        domain_aacode = []
        domain_num = []
        # protein code and respective number
        protein_aacode = []
        protein_num = []

        # fill d and p arrays with domain and protein sequences
        for i in range(0, delta):
            domain_aacode.append(l1[i])
            protein_aacode.append(l2[i])

        # compute the original positions in the system

        j1 = -1
        j2 = -1
        for i in range(0, len(domain_aacode)):
            if domain_aacode[i] != '.':
                j1 += 1
                domain_num.append(str(d_init + j1))
            if domain_aacode[i] == '.':
                domain_num.append('')
            if protein_aacode[i] != '-':
                j2 += 1
                protein_num.append(str(p_init + j2))
            if protein_aacode[i] == '-':
                protein_num.append('')
            output1.write(domain_aacode[i] + '\t' + str(domain_num[i]) + '\t' + protein_aacode[i] + '\t' + str(protein_num[i]) + '\n')

        output1.close()

    ## Matching
    pdb_2_pfam_dic = {}
    pfam_2_pdb_dic = {}

    for i in range(0, len(domain_num)):
        if domain_aacode[i] != '.' and protein_aacode[i] != '-':
            pdb_2_pfam_dic[int(protein_num[i])] = int(domain_num[i])
            pfam_2_pdb_dic[int(domain_num[i])] = int(protein_num[i])

    # outer loop runs along all DCA pairs
    N = len(dca_data)
    mapped_dca = []
    for i in range(N):
        try:
            dca_resi = dca_data['i'].iloc[i]
            dca_resj = dca_data['j'].iloc[i]
            di = dca_data['di'].iloc[i]
            diapc = dca_data['diapc'].iloc[i]
            mi = dca_data['mi'].iloc[i]
            if dca_resi and dca_resj in pfam_2_pdb_dic:
                mapped_dca.append([int(pfam_2_pdb_dic[dca_resi]), int(pfam_2_pdb_dic[dca_resj]), di, diapc, mi])
        except:
            print("%s" % dca_resi, dca_resj)
    x = np.array(mapped_dca)
    df_header = ["i", "j", "di", "diapc", "mi"]
    df_map = pd.DataFrame(x, columns=df_header)
    return df_map.astype({"i": "int", "j": "int", "di": "float", "diapc": "float", "mi": "float"})
