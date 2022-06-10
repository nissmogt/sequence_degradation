import numpy as np
from Bio import pairwise2
from mapping_functions import map_indices, apply_map

"""
Converts PDB indices.
Needs:
    - PDB sequences
    - PDB distance dataframe
    - PDB pairs to convert
Functions:
    - n = [to convert]
    - l1 = [pdbA pairs]
    - l2 = [pdbB pairs]
    - d = align(seqA, seqB)
    - m = map(l1, l2, d)
    - convert(n, m)
"""


def file_processing(file_in):
    labels = ["h1", "c1", "h2", "c2"]
    with open(file_in, "r") as f:
        header_1 = f.readline().rstrip()
        chain_1 = f.readline().rstrip()
        header_2 = f.readline().rstrip()
        chain_2 = f.readline().rstrip()

    return dict(zip(labels, [header_1, chain_1, header_2, chain_2]))


def align(sequence_1, sequence_2):
    # align(seq_a, seq_b)
    alignments = pairwise2.align.globalxs(sequence_2, sequence_1, -.5, -.1)
    print(pairwise2.format_alignment(*alignments[0]))
    map_reference = map_indices(alignments[0][0], 1, 0, alignments[0][1], 1, 0)
    # map_reference = map_reference.rename(columns={"i": "{}_i".format(pdb_name2), "A_i": "{}_res".format(pdb_name2),
    #                                               "j": "{}_i".format(pdb_name0), "A_j": "{}_res".format(pdb_name0)})
    return map_reference


# def run():
# File names
# f0 = "1R8U.fasta"
# f1 = "1L8C.fasta"
# f2 = "7LVS.fasta"
# pdb_name0 = f0.strip(".fasta")
# pdb_name1 = f1.strip(".fasta")
# pdb_name2 = f2.strip(".fasta")
# ff0 = "{}_FN_APC_CITED2_TAZ1_hmmer_inter_mapped_ca_dist.txt".format(pdb_name0)
# ff1 = "{}_FN_APC_TAZ1_HIF1A_hmmer_inter_mapped_ca_dist.txt".format(pdb_name1)
# map_reference_outfile = "reference_map_{}_{}_{}.txt".format(pdb_name2, pdb_name1, pdb_name0)
# map_reference_outfile = "reference_map_{}_{}.txt".format(pdb_name2, pdb_name0)
# pa = file_processing(f0)
# pb = file_processing(f1)
# pc = file_processing(f2)
#
# seq_a = pa["c1"] + pa["c2"]
# seq_b = pb["c1"]
# seq_c = pc["c1"] + pc["c2"][:27]
# Align sequences and save to file
# map_ref_hif = align(seq_b, seq_c)
# map_ref_cited = align(seq_a, seq_c)
# np.savetxt(map_reference_outfile, map_ref_cited, header="pdb_i\tpdb_res\tdca_i\tdca_res", fmt="%s\t%s\t%s\t%s",
#            comments='')
# return map_ref

# Make dictionary of mapping reference
# map_ref = run()
# map_pdb_dca_list = map_ref_cited.dropna()
# map_dict = dict(zip(map_pdb_dca_list["{}_i".format(pdb_name0)], map_pdb_dca_list["{}_i".format(pdb_name2)]))
#
# resi, resj = np.loadtxt(f0, unpack=True, usecols=(0, 1), skiprows=1, dtype=int, delimiter="\t")
# top_pairs = resi[:10]
# reindexed_pairs = np.zeros_like(top_pairs)
# for i, p in enumerate(top_pairs):
#     reindexed_pairs[i] = map_dict[str(p)]
