#!/home/kmm5/anaconda2/bin/python
"""
Author: Kareem Mehrabiani
Date: 17 October 2017

usage: calc_coverage.py [-h] [-gf] msa_file gap_threshold

Filters sequences of a Multiple-Sequence Alignment, given fraction of gaps (-)
present.

positional arguments:
  msa_file       Multiple-Sequence Alignment file
  gap_threshold  Gap threshold (from 0 to 1.0) Default=0.25

optional arguments:
  -h, --help     show this help message and exit
  -gf            Outputs a file containing the fraction of gaps for each
                 sequence.
"""

import argparse
import pandas as pd
import numpy as np
import re

outfile_ref = 'mapped_reference.txt'
outfile_map = 'mapped_match.txt'


def process_scan(scan_in):
    """Separates hmmscan results."""
    f = open(scan_in, 'r')
    data = f.read()
    raw_data = (re.findall(r'== domain (.*?)Internal', data, re.DOTALL))
    len(raw_data)

    scan = raw_data[0].split('\n')
    domain = scan[1].split()
    domain_seq = domain[2]
    pdb = scan[3].split()
    pdb_seq = pdb[2]
    return domain, pdb, domain_seq, pdb_seq


###############################################################################################
def map_loop(scan_in):
    """Creates a dictionary object with DCA-to-PDB map."""
    # Initialize
    file_header = 'Domain\td_seq\tp_seq\tPDB'
    domain, pdb, domain_seq, pdb_seq = process_scan(scan_in)
    domain_begin = int(domain[1])
    pdb_begin = int(pdb[1])
    len_domain = len(domain_seq)
    len_pdb = len(pdb_seq)
    domain_num = np.zeros(len_domain, dtype=int)
    pdb_num = np.zeros(len_pdb, dtype=int)
    resid_domain = []
    resid_pdb = []
    insert_count = 0
    gap_count = 0

    # Loop through domain and pdb sequence, counting every non-insert
    for i, res in enumerate(domain_seq, domain_begin):
        resid_domain.append(res)
        if res == '.':
            insert_count = insert_count + 1
        if res != '.':
            domain_num[i - int(domain[1])] = i - insert_count
    for i, res in enumerate(pdb_seq, pdb_begin):
        resid_pdb.append(res)
        if i > int(pdb[3]):
            if res == '-':
                gap_count = gap_count + 1
            if res != '-':
                pdb_num[i - int(pdb[1])] = (i + 4) - gap_count
        else:
            if res == '-':
                gap_count = gap_count + 1
            if res != '-':
                pdb_num[i - int(pdb[1])] = i - gap_count

    map_array = np.transpose([domain_num, resid_domain, resid_pdb, pdb_num])
    np.savetxt(outfile_ref, map_array, fmt='%s', delimiter='\t', header=file_header)

    map_dict = dict(zip(domain_num[domain_num != 0] - (domain_begin - 1),
                        pdb_num[pdb_num != 0] - (domain_begin - 1)))
    return map_dict


###############################################################################################
def map_dca(scan_in, dca_data):
    """Converts DCA pairs into PDB coordinates and outputs a file."""
    map_dict = map_loop(scan_in)
    N = len(dca_data)
    mapped_dca = []
    for i in range(N):
        try:
            dca_resi = dca_data['i'].iloc[i]
            dca_resj = dca_data['j'].iloc[i]
            di = dca_data['di'].iloc[i]
            diapc = dca_data['diapc'].iloc[i]
            mi = dca_data['mi'].iloc[i]
            if dca_resi and dca_resj in map_dict:
                # mapped_dca.append([int(map_dict[dca_resi]), int(map_dict[dca_resj]), di, diapc, mi,
                #                    int(dca_resi), int(dca_resj)])
                mapped_dca.append([int(map_dict[dca_resi]), int(map_dict[dca_resj]), di, diapc, mi])
        except:
            print("%s" % dca_resi, dca_resj)
    x = np.array(mapped_dca)
    df_header = ["i", "j", "di", "diapc", "mi"]
    # df_header = ["i", "j", "di", "diapc", "mi", "si", "sj"]
    df_map = pd.DataFrame(x, columns=df_header)
    return df_map.astype({"i": "int", "j": "int", "di": "float", "diapc": "float", "mi": "float"})
    # return df_map.astype(
    #     {"i": "int", "j": "int", "di": "float", "diapc": "float", "mi": "float", "si": "int", "sj": "int"})


###############################################################################################
def main():
    """
    Main function that calls .
    """
    import pandas as pd

    # Adds a description for the help section
    parser = argparse.ArgumentParser(description='')

    # Creates arguments for program
    parser.add_argument('scan_file', help='hmmscan output')
    parser.add_argument('dca_file', help='Three-column DCA pair list')
    parser.add_argument('--sep', type=int, default=4,
                        help='Sequence separation (default: 4)')
    args = parser.parse_args()

    # Input arguments for gap filter function
    scan_in = args.scan_file
    dca_in = args.dca_file
    separation = args.sep
    print("\nSequence sep: %d\n" % separation)

    names = ['residue_i', 'residue_j', 'scores']
    dca_data = pd.read_csv(dca_in, names=names, delimiter=',', usecols=(0, 1, 2))
    dca_data['dist'] = np.abs(dca_data['residue_i'].sub(dca_data['residue_j'], axis=0))
    # Keep sequences with |i-j|>=separation
    dca_data = dca_data[dca_data['dist'] > separation]
    # Rank pairs by DCA score
    dca_data = dca_data.sort_values(by=['scores'], ascending=False)

    ### MAP DCA ###
# fin = "D:\\bootstrap\\PF00022\\replicates\\sub0\\mf\\pc0.2\\DI_PF00022_n19653.txt"
# dca_data = pd.read_csv(fin, delimiter=',')
# df=map_dca("2zwh_scan.txt", dca_data)

# if __name__ == '__main__':
#     main()
