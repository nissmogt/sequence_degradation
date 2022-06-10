import numpy as np

msa = "TAZ1_CITED2_hmmer.fas"
# msa = "cited2_hmmer_1r8u_filtered_25.fasta"
# msa = "taz1_7lvs_filtered_25.fasta"
entry = []
with open(msa, 'r') as fp:
    for i, line in enumerate(fp.readlines()):
        entry.append(line.rstrip())

headers = entry[0:len(entry):2]
sequences = entry[1:len(entry):2]

ai_count = 0
aj_count = 0
aij_count = 0
ai = 38
aj = 144
Ai = 'N'
Aj = 'Q'
seq_index = []
for k in range(len(sequences)):
    if sequences[k][ai] == Ai:
        ai_count += 1
    if sequences[k][aj] == Aj:
        aj_count += 1
    if sequences[k][ai] == Ai and sequences[k][aj] == Aj:
        aij_count += 1
        seq_index.append(k)
# aij_count = ai_count
for j in range(len(seq_index)):
    print("SEQ:{} HEAD:{}  COUNT:{}".format(seq_index[j], headers[seq_index[j]], aij_count))
