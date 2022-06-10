import numpy as np
import pandas as pd


def pipeline_pdb_distance_matrix(pdbid, cutoff, heavy_atom=True, read=False, plot=False):
    """
    Used in calculating distance matrix for given msa input.
    :param pdbid: string - Name for PDB file should follow the following format 'PDBID_CHAIN1_PDBID_CHAIN2'
    :param heavy_atom: string - 'aa' for all atom min distance calculation. 'ca' for c-alpha distance calculation.
    :param cutoff: float - Cutoff defines contact between two inter-chain residues in PDB.
                   Note: 8A is hard-coded for monomer.
    :param read: boolean - If True, reads distance matrix file into pandas Dataframe. (Default: False)
    :param plot: boolean - Plot distance matrix at chosen cutoff using chosen cutoff type. (Default: False)
    :return: List of two pandas Dataframes. Monomer pairs and interface pairs.
    """
    from distance_pdb import distance_matrix
    from import_distance_matrix import import_pdb_distance_matrix
    from plot_contact_map import plot_cm

    print("\n\t-- |{}| DISTANCE MATRIX CALCULATION at |{}| inter-cutoff: |{}| --".format(pdbid, heavy_atom, cutoff))
    if read:
        df_pdb = import_pdb_distance_matrix(pdbid, heavy_atom=heavy_atom)
    else:
        df_pdb, ch = distance_matrix(pdbid, heavy_atom=heavy_atom)

    chain1 = np.unique(df_pdb.chain_1.values)[0]
    chain2 = np.unique(df_pdb.chain_1.values)[1]
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    total_length = max(df_mon[df_mon["chain_2"] == chain2].j)
    chain1_length = max(df_mon[df_mon["chain_1"] == chain1].j)
    chain2_length = total_length - chain1_length

    if heavy_atom:    # monomer definition depends on atom-atom resolution
        mt = 5.0
    else:
        mt = 8.0
    df_mon = df_mon[df_mon["d"] <= mt]    # hard-coded monomer cutoff
    df_inter = df_inter[df_inter["d"] <= cutoff]

    # total_length = sum(chain_lengths)
    print("\t||Chain 1\t||Chain 2\nlengths: {}\t||{}\t||Total length: {}".format(chain1_length, chain2_length,
                                                                                 total_length))
    pdb_df_list = [df_pdb, df_mon, df_inter]
    if plot:
        df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm
        plot_cm(pdb_df_list[1:], cutoff=cutoff, length_a=chain1_length, length=total_length, atom=heavy_atom,
                df_dca=df_empty, msa_name=pdbid)

    return pdb_df_list, [chain1_length, chain2_length, total_length]


def pipeline_mapping(msa_name, df_dca, flag, read=False, mfdca=False):
    """
    Map DCA indices to PDB-distance-matrix indices
    :param mfdca:
    :param flag: Boolean, set to True if using DBMarks a2m files
    :param read: Boolean, either read from a reference file or generate reference file (default).
    :param df_dca: Pandas Dataframe, DCA input
    :param msa_name: String, Name of msa in ID_CHAIN1_ID_CHAIN2 format
    :return: Pandas Dataframe, mapped DCA output
    """
    from msa_functions import read_target_sequence_in_msa
    from read_db import get_lengths
    from get_residues import get_residues
    from get_region import get_dca_indices
    from mapping_functions import align_dca2pdb, apply_map, apply_map_mfdca
    print("(pipeline mapping)")
    if read:
        infile = "ref_map_{}.txt".format(msa_name.strip(".fas"))
        map_pdb_dca = pd.read_csv(infile, delimiter="\t", header=0, dtype=str)
        # map_pdb_dca = map_pdb_dca.replace("?", np.nan).dropna()    # some pdbs have unknown seq res UNK
        # map_pdb_dca = pd.read_csv(infile, delimiter="\t", header=0, dtype=str, usecols=(0,1)) # used for HK-RR
        # map_pdb_dca["#HMM"] = map_pdb_dca["#HMM"].astype(int) + 1 # used for HK-RR
        # map_to_pdb = dict(zip(map_pdb_dca["#HMM"], map_pdb_dca["col"])) # used for HK-RR
        map_pdb_dca = map_pdb_dca.replace("X", np.nan).dropna()    # some pdbs have unknown seq res UNK
        map_to_pdb = dict(zip(map_pdb_dca["dca_i"], map_pdb_dca["pdb_i"]))

    else:
        uniprot_lengths = get_lengths(msa_name)
        if flag:
            _, dca_lengths, _ = get_dca_indices(msa_name, uniprot_lengths[0])
        else:
            dca_lengths = uniprot_lengths
        # -- GET MAP FROM MSA TO PDB --
        pdbseq_1, pdbseq_2 = get_residues(msa_name, seq=True)
        pdbseq = [pdbseq_1, pdbseq_2]
        # splits msa sequence based on modified uniprot lengths (removed lowercase)
        msaseq = read_target_sequence_in_msa(msa_name, split=True, len_a=dca_lengths[0])
        map_to_pdb = align_dca2pdb(msa_name, pdbseq, msaseq)

    # print("(map dictionary) {}".format(map_to_pdb))
    if mfdca:
        mapped_dca_array = apply_map_mfdca(df_dca.to_numpy(), map_to_pdb)
    else:
        mapped_dca_array = apply_map(df_dca.to_numpy(), map_to_pdb)
    return mapped_dca_array


def pipeline_interface(df_dca_mapped, pdb_chain1_length):
    """

    :param df_dca_mapped: 
    :param pdb_chain1_length:
    :return:
    """
    df_dca_mapped_inter = df_dca_mapped.loc[(df_dca_mapped["i"] < pdb_chain1_length) &
                                            (df_dca_mapped["j"] > pdb_chain1_length)]
    return df_dca_mapped_inter
