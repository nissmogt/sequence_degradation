# pdb = "structures\\7lvs.pdb"
# all_atom = True
def distance_matrix(pdb_in, heavy_atom=False):
    """
    Calculates distance matrix and outputs it to a csv file.
    :param pdb_in:String; Full PDB file path and name.
    :param heavy_atom: Boolean; Calculate heavy atom distance or Calpha.
    :return: Dataframe and List of chain lengths
    """
    import time
    import os
    import pandas as pd
    from itertools import combinations_with_replacement
    from get_residues import get_residues
    from distance_functions import calc_min_dist, calc_ca_distance
    from parse_pdb import parse_pdb

    # Checks if pdb_in has a directory included, if not set output dir to working dir
    if os.path.dirname(pdb_in):
        pdb = os.path.basename(pdb_in)
        pdb_dir = os.path.dirname(pdb_in)
    else:
        pdb = pdb_in
        pdb_dir = ""
    # array initialization
    resi_list = []
    resj_list = []
    actual_i_list = []
    actual_j_list = []
    distance_list = []
    chain_1_list = []
    chain_2_list = []
    atom_id_list = []
    residue_list = []

    # Function information output
    print("\t(distance matrix) begin loop")
    chain_model = parse_pdb(pdb_in)

    residues, chain_lengths = get_residues(chain_model)  # Output list of residues from pdb
    # count every pair of residues NOTE: INDEX BEGINS AT 0 BUT 1 IS ADDED BELOW
    pair_list = combinations_with_replacement(range(len(residues)), 2)

    start_time = time.time()
    for i, j in pair_list:
        if i != j:  # ensure residue i not equal to j
            res_a = residues[int(i)]
            res_b = residues[int(j)]
            actual_i_list.append(res_a.id[1])
            actual_j_list.append(res_b.id[1])
            # get chain id
            chain_1_list.append(res_a.get_parent().id)
            chain_2_list.append(res_b.get_parent().id)
            # resets res index to 1 SEE NOTE ABOVE.
            resi_list.append(i + 1)
            resj_list.append(j + 1)
            residue_list.append((res_a.resname, res_b.resname))
            if heavy_atom:
                mindist, atom_ids = calc_min_dist(res_a, res_b)
                distance_list.append(mindist)
                atom_id_list.append(atom_ids)
            else:
                if res_a.has_id("CA") and res_b.has_id("CA"):
                    distance_list.append(calc_ca_distance(res_a, res_b))
                else:
                    print(f"NOTE:Res {res_a.get_full_id()}\n\tor {res_b.get_full_id()} not calculated! (missing CA)\n")
    # fileout.close()
    print("\t -- LOOP TIME -- {}".format(time.time() - start_time))
    # makes a pandas dataframe
    if heavy_atom:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list, 'chain_1': chain_1_list,
                               'chain_2': chain_2_list, 'resnames': residue_list, 'atom_id': atom_id_list})
        filename = os.path.join(pdb_dir, f"atom_distance_matrix_{pdb.split('.pdb')[0]}.txt")
        header = "i\tj\tdist_aa\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id"
        df_pdb.to_csv(filename, sep='\t', index=False, header=header, float_format='%.5f')
    else:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list,
                               'chain_1': chain_1_list, 'chain_2': chain_2_list, 'resnames': residue_list})
        filename = os.path.join(pdb_dir, f"ca_distance_matrix_{pdb.split('.pdb')[0]}.txt")
        header = "i\tj\tdist_ca\tsi\tsj\tchain_1\tchain_2\tresnames"
        df_pdb.to_csv(filename, sep='\t', index=False, header=header, float_format='%.5f')

    print(f"wrote {filename}")
    return df_pdb, chain_lengths
