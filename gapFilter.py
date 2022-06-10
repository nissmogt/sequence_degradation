import os.path

import numpy as np
import csv


def gap_filter(msa_input, gap_threshold=0.25):
    """
    First removes insertions (.) and lowercase letters, if there are any.
    Then counts the number of gaps in each sequence and from this,
    calculates the fraction of gaps present. If this fraction is greater
    than the input gap threshold (Default=0.25), then that
    particular sequence is thrown out. Optional: Gap percent per sequence
    can be returned (by output_gaps=True).

    Author: Kareem Mehrabiani
    Date: 17 October 2017
    Updated 21 October 2020
    """

    from Bio import AlignIO
    from Bio import SeqIO
    from Bio.Seq import Seq
    print('================================================================')
    print('Use gapFilter.py -h or --help for usage.')
    print('Gap threshold is {}%.'.format(gap_threshold * 100))
    q = msa_input.split("\\")
    if len(q) > 1:
        msaDir = os.path.dirname(msa_input) + "\\"
        msa_input = os.path.basename(msa_input)
        print(f"CHECK OK. DIR:{msaDir}\tPFAM:{msa_input}")
    else:
        msaDir = ""  # change this to your directory

    alignment = AlignIO.read('{}{}'.format(msaDir, msa_input), 'fasta')
    total_sequences = len(alignment)

    output_handle = open("{}{}_filtered_{}.fasta".format(msaDir, msa_input.strip(".txt"),
                                                         int(gap_threshold * 100)), 'w', encoding='utf-8')
    # removed_output = open('removed_seqid' + '_filtered_%dp.txt' % (gap_threshold * 100), 'w')
    print('Number of sequences in MSA: %d\n' % total_sequences)

    filtered_sequences = 0
    percent_gaps = np.zeros(total_sequences)

    for seq_num, record in enumerate(alignment):
        # Counts number of gaps in current sequence
        # Note: removes insertions via ungap method
        current_sequence = record.seq.ungap('.')
        temp_seq = Seq('')

        # Removes lowercase letters
        for residue in range(len(current_sequence)):
            if not current_sequence[residue].islower():
                temp_seq += current_sequence[residue]

        current_sequence = temp_seq

        # Calculates % gaps in current sequence
        gap_count = current_sequence.count('-')
        percent_gaps[seq_num] = float(gap_count) / float(len(current_sequence))

        # Outputs current sequence that maintains minimum gap threshold
        if percent_gaps[seq_num] <= gap_threshold:
            filtered_sequences += 1
            output_handle.write('>%s\n' % record.id + ''.join(current_sequence) + '\n')
        # else:
        #     removed_output.write('%s\n' % record.id)

    removed_sequences = total_sequences - filtered_sequences

    print('Number of sequences removed: %d(%s%%)' % (
        removed_sequences, removed_sequences * 100 / total_sequences))
    print('Number of sequences kept: %d\n' % filtered_sequences)
    print('Wrote new file: \'' + output_handle.name + '\'\n')
    output_handle.close()
    # removed_output.close()

    return output_handle.name, filtered_sequences, percent_gaps


# msaName = '1EM8_D_1EM8_C.fas'
def batch_filter(list_msa):
    gaps = 0.25
    nTotal = []
    with open("benchmark_systems_FNi_count182.txt", "r", encoding='utf-8') as msaList:
        for msa in msaList:
            msa_name = msa.rstrip("\n")

            # Gap filter function call ###
            nKept, percentGaps = gap_filter(msa_name, gaps)
            nTotal.append([msa_name, nKept])

            # Write gap fractions per sequence to a file
            gapFile = 'gap_fraction\\{}_gap_fraction_{}.txt'.format(msa_name.strip(".fas"), int(gaps * 100))
            np.savetxt(gapFile, percentGaps, fmt='%.3f', newline='\n')
            print('Also wrote \'' + gapFile + '\'\n')

        outfile = "filtered_25p_{}.csv".format(msaList.name.strip(".txt"))
        with open(outfile, "w", encoding='utf-8', newline='') as output:
            writer = csv.writer(output)
            writer.writerow(["MSA", "nSeq_filtered_25"])
            writer.writerows(nTotal)

# msa_list = "benchmark_systems_FNi_count182.txt"
# batch_filter(msa_list)
