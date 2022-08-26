import os
import sys
import pandas as pd


def parse_pdb(pdb_in):
    import Bio.PDB
    if os.path.dirname(pdb_in):
        dir_pdb = os.path.dirname(pdb_in)
        pdb_file = os.path.basename(pdb_in)
    else:
        dir_pdb = ""
        pdb_file = pdb_in

    pdb_id, ext = os.path.splitext(pdb_file)
    chain_lengths = []
    parser = Bio.PDB.PDBParser()
    struct = parser.get_structure(pdb_id, pdb_in)
    model = struct[0]
    chains = [ch_id for ch_id in model.get_list()]
    return chains


def calc_min_dist(res_a, res_b):
    import numpy as np
    na_atoms = len(res_a.get_list())
    nb_atoms = len(res_b.get_list())
    atomic_matrix = np.zeros((na_atoms, nb_atoms))
    for i, atom_i in enumerate(res_a):
        for j, atom_j in enumerate(res_b):
            if atom_i.mass > 2 and atom_j.mass > 2:
                atomic_matrix[i][j] = np.linalg.norm(atom_i.get_coord() - atom_j.get_coord())
    mindist = np.amin(atomic_matrix[np.nonzero(atomic_matrix)])
    atom_i_index, atom_j_index = np.where(atomic_matrix == mindist)  # gives index of mindist atoms i,j
    min_atom_i = res_a.get_list()[atom_i_index[0]]
    min_atom_j = res_b.get_list()[atom_j_index[0]]
    return mindist, (min_atom_i.id, min_atom_j.id)


def calc_ca_distance(res_a, res_b):
    """
    Calculates the distance between a pair of CA atoms
    :param res_a: Biopython residue object - residue a
    :param res_b: Biopython residue object - residue b
    :return: Distance between CA atoms
    """
    import numpy as np
    a = res_a["CA"].get_coord()
    b = res_b["CA"].get_coord()
    dist = np.linalg.norm(a - b)
    return dist


def three2one(sequence):
    """ Lookup table - translate a protein sequence from 3 to 1 letter code
    """

    code = {"GLY": "G", "ALA": "A", "LEU": "L", "ILE": "I",
            "ARG": "R", "LYS": "K", "MET": "M", "CYS": "C",
            "TYR": "Y", "THR": "T", "PRO": "P", "SER": "S",
            "TRP": "W", "ASP": "D", "GLU": "E", "ASN": "N",
            "GLN": "Q", "PHE": "F", "HIS": "H", "VAL": "V",
            "M3L": "K", "MSE": "M", "CAS": "C"}

    newprot = ""
    for a in sequence:
        newprot += code.get(a, "?")

    return newprot


def get_residues(chains, seq=False):
    """
    Build a simple list of residues from a single chain of a PDB file.
    :param chains: Bio PDB Chain Object
    :param seq: Boolean (Default: False) - Outputs sequence if True.
    :return: A list of Bio.PDB.Residue objects.
    """
    import Bio.PDB
    # from parse_pdb import parse_pdb
    # chains = parse_pdb(pdb)
    chain_lengths = []
    residues = []
    sequence = []
    for ch in chains:
        # make sure res are standard AA
        num_residues = 0
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            # if Bio.PDB.is_aa(res, standard=True):
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            num_residues += 1
            if (
                    res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS' or res_id == 'HMS') and is_regular_res:
                residues.append(res)
                sequence.append(res.get_resname())
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))
        chain_lengths.append(num_residues)

    if seq:
        sequence = three2one(sequence)
        seq_a = sequence[:chain_lengths[0]]
        seq_b = sequence[chain_lengths[0]:]
        return seq_a, seq_b
    else:
        return residues, chain_lengths


def get_residues_web(pdb, seq=False, chain_ids=None):
    """
    Build a simple list of residues from a single chain of a PDB file.
    :param pdb: String - PDB filename
    :param seq: Boolean (Default: False) - Outputs sequence if True.
    :param chain_ids:
    :return: A list of Bio.PDB.Residue objects.
    """
    import Bio.PDB

    if chain_ids is None:
        chain_ids = ["A", "B"]

    # chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]
    chain_lengths = []
    parser = Bio.PDB.PDBParser()

    struct = parser.get_structure(pdb, pdb)
    model = struct[0]
    # if len(self.chain_ids) == 0:
    # get residues from every chain.
    #    chains = model.get_list()
    # else:
    chains = [model[ch_id] for ch_id in chain_ids]

    residues = []
    sequence = []
    for ch in chains:
        # make sure res are standard AA
        num_residues = 0
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            # if Bio.PDB.is_aa(res, standard=True):
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
                residues.append(res)
                sequence.append(res.get_resname())
                num_residues += 1
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))
        chain_lengths.append(num_residues)

    if seq:
        sequence = three2one(sequence)
        seq_a = sequence[:chain_lengths[0]]
        seq_b = sequence[chain_lengths[0]:]
        return seq_a, seq_b
    else:
        return residues, chain_lengths


def distance_matrix(pdb_in, heavy_atom=False):
    """
    Calculates distance matrix and outputs it to a csv file.
    :param pdb_in:String; Full PDB file path and name.
    :param heavy_atom: Boolean; Calculate heavy atom distance or Calpha.
    :return: Dataframe and List of chain lengths
    """
    import time
    from itertools import combinations_with_replacement

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
    dir_matrix = os.path.join(pdb_dir, "matrix")
    if not os.path.exists(dir_matrix):
        os.makedirs(dir_matrix)
    if heavy_atom:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list, 'chain_1': chain_1_list,
                               'chain_2': chain_2_list, 'resnames': residue_list, 'atom_id': atom_id_list})
        filename = os.path.join(dir_matrix, f"atom_distance_matrix_{pdb.split('.pdb')[0]}.txt")
        header = "i\tj\tdist_aa\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id"
        df_pdb.to_csv(filename, sep='\t', index=False, header=header, float_format='%.5f')
    else:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list,
                               'chain_1': chain_1_list, 'chain_2': chain_2_list, 'resnames': residue_list})
        filename = os.path.join(dir_matrix, f"ca_distance_matrix_{pdb.split('.pdb')[0]}.txt")
        header = "i\tj\tdist_ca\tsi\tsj\tchain_1\tchain_2\tresnames"
        df_pdb.to_csv(filename, sep='\t', index=False, header=header, float_format='%.5f')

    print(f"wrote {filename}")
    return df_pdb, chain_lengths


def pipeline_pdb(pdb_id, dir_pdb):
    # PDB distance matrix
    # PDB_id = "5pti"
    # PDB_id = "1or7"
    matrix_file = os.path.join(dir_pdb, "matrix", f"atom_distance_matrix_{pdb_id}.txt")
    if os.path.exists(matrix_file):
        print("Path exists! Reading from file instead...")
        pdb_dataframe = pd.read_csv(matrix_file, header=0, delimiter="\t")
    else:
        # uncomment if not reading
        pdb_path = os.path.join(dir_pdb, f"{pdb_id}.pdb")
        pdb_dataframe, chain_len = distance_matrix(pdb_path, heavy_atom=True)

    # min_i = pdb_dataframe[pdb_dataframe.si > 0]["i"].min()
    # pdb_dataframe["i"] = pdb_dataframe["i"] - (min_i - 1)
    # pdb_dataframe["j"] = pdb_dataframe["j"] - (min_i - 1)
    # out_df = pdb_dataframe[pdb_dataframe.si > 0]
    # return out_df.reset_index(drop=True)
    return pdb_dataframe
