import time
import os
import numpy as np
from Bio import Align, AlignIO


def test_suite():
    test_seq = np.array([[1, 2, 2, 10], [1, 2, 2, 10], [1, 1, 3, 10]])
    assert effective_seq_calc(test_seq) == 2
    return "All Good."


def aa2num(aa):
    """convert aa into num"""
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    states = len(alphabet)
    a2n = {}
    for a, n in zip(alphabet, range(states)):
        a2n[a] = n
    if aa in a2n:
        return a2n[aa]
    else:
        return a2n['-']


def letter2number(a):
    d = {'-':
             1,
         'A':
             2,
         'C':
             3,
         'D':
             4,
         'E':
             5,
         'F':
             6,
         'G':
             7,
         'H':
             8,
         'I':
             9,
         'K':
             10,
         'L':
             11,
         'M':
             12,
         'N':
             13,
         'P':
             14,
         'Q':
             15,
         'R':
             16,
         'S':
             17,
         'T':
             18,
         'V':
             19,
         'W':
             20,
         'Y':
             21,
         'X':
             1,
         'U':
             1,
         'O':
             1
         }
    if a not in d:
        return 1
    else:
        return d[a]


def read_target_sequence_in_msa(msa_name, split=False, len_a=None):
    """
    Reads first sequence in MSA which should be a template
    :return: List - MSA sequence
    """
    # Create MSA-seq-template file object
    msa_dir = ""
    msa_file_obj = open("{}{}.fasta".format(msa_dir, msa_name), 'r')
    msa_header = msa_file_obj.readline().rstrip()
    msa_seq = msa_file_obj.readline().rstrip()
    msa_file_obj.close()
    if split:
        assert len_a < len(msa_seq)
        return msa_seq[:len_a], msa_seq[len_a:]
    else:
        return msa_seq


def parse_fasta(msa_name, skip_first=True, limit=-1):
    """
    function to parse fasta
    :return header and sequence of fasta
    """
    import os
    header = []
    sequence = []
    # msa_path = "PDB_benchmark_alignments\\dimers\\"
    msa_path = "../"

    if os.path.exists(msa_path):

        # lines = open("{}{}".format(msa_path, msa_name), 'r')
        lines = open(msa_name, 'r')
        if skip_first:
            # used to skip first fasta sequence and header
            skip_null = [next(lines) for x in range(2)]
        for line in lines:
            line = line.rstrip()
            if line[0] == ">":
                if len(header) == limit:
                    break
                header.append(line[1:])
                sequence.append([])
            else:
                sequence[-1].append(line)
        lines.close()
        sequence = [''.join(seq) for seq in sequence]
        return np.array(header), np.array(sequence)


def mon_neff(sys_name):
    from scipy.spatial.distance import pdist, squareform
    from scramble_sequence import split_header_seq
    from read_db import get_lengths
    from get_region import get_dca_indices
    from aa2num import aa2num

    ch = get_lengths(sys_name)
    _, dca_lengths, _ = get_dca_indices(sys_name, ch[0])
    h1, h2, s1, s2, s = split_header_seq(sys_name, dca_lengths[0], full=True)
    ss = [s1, s2, s]

    neff = []
    for si in ss:
        print(si)
        msa_in_number_form = []
        for seq in si:
            msa_in_number_form.append([aa2num(aa) for aa in seq])
        msa_in_number_form = np.array(msa_in_number_form)

        x = pdist(msa_in_number_form, metric='hamming')
        y = 1 - squareform(x)
        w = (y >= 0.8).astype(np.float)
        w = 1 / np.sum(w, -1)
        neff.append(np.sum(w))


def effective_seq_calc(sequence, eff_cutoff=0.8):
    """compute effective weight for each sequence"""
    from scipy.spatial.distance import pdist, squareform

    # pairwise identity
    msa_sm = 1.0 - squareform(pdist(sequence, "hamming"))
    # weights = 1. / (1 + np.sum(squareform(pdist(Y, "hamming")), axis=0))

    # weight for each sequence
    msa_w = (msa_sm >= eff_cutoff).astype('float')
    msa_w = 1 / np.sum(msa_w, -1)
    return np.sum(msa_w)


def remove_randseq(x, seed):
    n_seqs = len(x)
    rng = np.random.default_rng(seed)
    n_to_keep = n_seqs // 2
    to_keep = rng.choice(n_seqs, size=n_to_keep, replace=False)
    # new_msa = x[idx]
    # ratio = len(idx) / n_seqs
    # print("YO {}".format(ratio))
    return to_keep


def msa_object(msaName, skip_first=False, all_rw=False):
    """converts list of sequences to msa"""

    print("parsing msa...")
    headers, sequences = parse_fasta(msaName, skip_first=skip_first)

    msa_in_number_form = []
    print("encoding to number form...")
    for ii, seq in enumerate(sequences):
        msa_in_number_form.append([letter2number(aa) for aa in seq])
    msa_in_number_form = np.array(msa_in_number_form)

    # remove positions with more than > 50% gaps
    # msa, v_idx = self.filt_gaps(msa_ori, 0.5)
    msa = msa_in_number_form

    return msa, [headers, sequences]


def msa_info(msa, all_rw=False):
    # compute effective weight for each sequence
    print("calculating effective sequences...")
    if all_rw:
        neff = []
        for i in range(1, 9):
            print("rw:{}".format(i))
            msa_weights = effective_seq_calc(msa, i * 0.10)
            neff.append(np.sum(msa_weights))
        ncol = msa.shape[1]  # length of sequence
        return {"neff": neff, "nrow": msa.shape[0], "ncol": ncol}
    else:
        msa_weights = effective_seq_calc(msa, 0.8)
        neff = (np.sum(msa_weights))
        ncol = msa.shape[1]  # length of sequence
        return {"neff": neff, "nrow": msa.shape[0], "ncol": ncol}


def count_sequence(msaName):
    """
    Reads msa fasta file and counts number of '>' characters which equals the number of sequences
    :param msaName:
    :return: int Number of sequences in fasta file
    """
    msaDir = "PDB_benchmark_alignments\\"
    with open("{}{}".format(msaDir, msaName), "r", encoding='utf-8') as msa:
        count = msa.read().count(">")
    return count


def stats_for_msa_list(msa_list):
    """
    For every msa in msa list count number of sequences, calculate number of effective sequences,
    and add these values to a list with the msa name.
    :param msa_list:
    :return: list containing MSA name, number of sequences, number of effective sequences
    """
    import pandas as pd
    stats_list = []
    start = time.time()
    for msa in msa_list:
        msaName = msa + '.fas'
        print(msaName)
        pdb_data = "distance_matrix\\heavy_atom_distance_matrix_{}.txt".format(msa)
        df_pdb = pd.read_csv(pdb_data, delimiter="\t")
        df_pdb = df_pdb[df_pdb["d"] <= 12]
        dfpdb = df_pdb[df_pdb["chain_1"] != df_pdb["chain_2"]].reset_index(drop=True)
        n_contacts = len(dfpdb)
        msaObject = msa_object(msaName, skip_first=True)
        nSeqs = msaObject["nrow"]
        nEff = msaObject["neff"]
        seqLength = msaObject["ncol"]
        stats_list.append([msaName.strip('.fas'), nEff / seqLength, n_contacts / seqLength, nSeqs, nEff, seqLength])
    print("Duration of loop: {}".format(time.time() - start))
    return stats_list


def get_msa_len(msa_data):
    with open(msa_data, "r") as msafile:
        Lines = msafile.readlines()
        L = len(Lines[1].rstrip())
    msafile.close()
    return L


def write_msa_stats(msa_list):
    import csv
    stats_list = stats_for_msa_list(msa_list)
    outfile = "dimers_stats.csv"
    with open(outfile, "w", encoding='utf-8', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(["MSA", "Neff/L", "inter-contact density", "Nseqs", "Neff", "L"])
        writer.writerows(stats_list)
    return stats_list


# msafiles = "benchmark_dimer_systems.txt"
# s = write_msa_stats(msafiles)

def make_monomer_msa_from_dimer(msa):
    from scramble_sequence import split_header_seq
    from read_db import get_lengths
    from get_region import get_dca_indices
    outDir = "monomer_alignments\\"
    # msa = "2OXG_Z_2OXG_Y"
    cid = [msa.split("_")[1], msa.split("_")[3]]
    ch = get_lengths(msa)
    _, dca_chains, _ = get_dca_indices(msa, length_a=ch[0])
    x = split_header_seq(msa, dca_chains[0])
    for i in range(2):
        a = np.array([x[i], x[i + 2]])
        a = a.transpose()
        np.savetxt("{}{}_{}.fas".format(outDir, msa[:4], cid[i]), a, fmt=">%s\n%s")


# dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
#       "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
#       '5MU7_B_5MU7_A']
# dimers = ['1EFP_A_1EFP_B', '1EP3_A_1EP3_B', '1ZUN_A_1ZUN_B', '3A0R_A_3A0R_B',
#           '3MML_A_3MML_B', '1B70_A_1B70_B', '1BXR_A_1BXR_B']

# from make_dimer_list import make_dimer_list
# dimers = make_dimer_list()
# s = write_msa_stats(dimers)
# make_monomer_msa_from_dimer(msaname)
# x = msa_object("TAZ1_HIF1A_hmmer.fas", skip_first=False)
# y = msa_object("TAZ1_CITED2_hmmer.fas", skip_first=False)
# x, ll = msa_object("msa_len252.fasta")
# msa_weights = effective_seq_calc(x, 0.8)
# neff = (np.sum(msa_weights))


def write_msa(seq_entry, output_filepath):
    directory = os.path.dirname(output_filepath)
    filename = os.path.basename(output_filepath)
    os.makedirs(directory, exist_ok=True)
    with open(output_filepath, "w") as output_handle:
        AlignIO.write(Align.MultipleSeqAlignment(seq_entry), output_handle, "fasta")


def generate_replicates(_msa, n_subalign, _dir_sub_out):
    import time
    import os

    msaDir = os.path.dirname(_msa)
    msa = os.path.basename(_msa)
    msaname = msa.split(".")[0]
    print(f"CHECK OK. DIR:{msaDir}\tPFAM:{msa}")
    # dir_sub_out = os.path.join(msaDir, "replicates")
    pid = os.path.basename(_msa).split("_")[0]
    print(f"CHECK OK. DIR:{_dir_sub_out}\tFilename:{_msa}")

    # Extract pfam id from filename and use to build folder structure
    # DIR_SUB_OUT = f"bootstrap\\{pid}\\subalignments\\"
    if not os.path.exists(_dir_sub_out):
        os.makedirs(_dir_sub_out)

    # Read in msa and insert into a list of decreasing nseq
    msa = AlignIO.read(open(_msa, "r"), "fasta")
    records = [record for record in msa]
    seed_list = []
    max_nseq = 100000
    for i in range(n_subalign):
        seed_list.append(time.time_ns())
        row = 0
        msa_sets = [records]
        nseq = len(msa_sets[0])
        if nseq < max_nseq:
            list_of_len = [nseq]
        else:
            list_of_len = []
        while nseq > 0:
            # Remove random rows from MSA until nseq = 1
            current_msa = msa_sets[row]
            idx_to_keep = remove_randseq(current_msa, seed_list[i])
            new_msa = [current_msa[ind] for j, ind in enumerate(idx_to_keep)]
            dir_sub_n = os.path.join(_dir_sub_out, f"sub{i}")
            fasta_out = os.path.join(dir_sub_n, f"{pid}_n{nseq}_sub{i}.fasta")
            msa_sets.append(new_msa)
            if nseq < max_nseq:
                write_msa(msa_sets[row], fasta_out)

            print(f"Number of seqs: {nseq}")
            row += 1
            nseq = len(new_msa)
            if nseq < max_nseq:
                list_of_len.append(nseq)
    assert len(np.unique(seed_list, return_counts=True)[0]) == len(seed_list)
    seed_lout = os.path.join(_dir_sub_out, "seed_list.txt")
    with open(seed_lout, "w") as out_handle:
        for k in range(len(seed_list)):
            out_handle.write(f"seed{k}:{seed_list[k]}\n")

    len_lout = os.path.join(_dir_sub_out, "length_list.txt")
    with open(len_lout, "w") as out_handle:
        for k in range(len(list_of_len[:-1])):
            out_handle.write(f"{list_of_len[k]}\n")

    return list_of_len


# for i in range(len(msa_sets)):
#         y, ll = msa_object(fp)
#         neff = effective_seq_calc(y, 0.8)
#         print(neff)

def test():
    fin = "testid_test.fasta"
    generate_replicates(fin)
# fin = "PF00014_filtered_25.fasta"
# fin = "PF04542_filtered_25.fasta"

# for j in range(11):
#     msa_weights = effective_seq_calc(msa_sets[j+3], 0.8)
#     neff.append(np.sum(msa_weights))
#     print(len(msa_sets[j+3]), neff[j])

# x, s = msa_object("testid_test.fasta", skip_first=False)
# mw = effective_seq_calc(x)
