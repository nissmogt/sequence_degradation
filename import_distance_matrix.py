def import_pdb_distance_matrix(pdb_file, heavy_atom=False):
    """
    :param pdb_file:
    :param heavy_atom: Boolean {Default: False}, specify coarseness of contact map
    :return: Dataframe, pdb pairs and other information of pairs
    """
    import pandas as pd
    import os

    out_path = "contact_maps\\"
    fname = "(import distance matrix)"
    pdbid = os.path.basename(pdb_file).split(".")[0]
    if heavy_atom:
        filename = "{}atom_distance_matrix_{}.txt".format(out_path, pdbid)
    else:
        filename = "{}ca_distance_matrix_{}.txt".format(out_path, pdbid)

    df_pdb = pd.read_csv(filename, delimiter="\t")
    return df_pdb


