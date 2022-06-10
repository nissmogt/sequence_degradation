import pandas as pd
import matplotlib.pyplot as plt
from rank_hamming import rank_hamming

pf = "taz1_avg_distance_matrix_8A.txt"
dfpdb = pd.read_csv(pf, delim_whitespace=True)
f = "FN_APC_PF02135_TAZ1_filtered_25.fasta.txt"
df = pd.read_csv(f, delimiter=',', names=['i', 'j', 'fn', 'cn'])
df_ranked = rank_hamming(df, 'cn')
n = 10
plt.scatter('i', 'j', data=dfpdb)
plt.scatter('i', 'j', data=df_ranked[:n])
plt.show()
