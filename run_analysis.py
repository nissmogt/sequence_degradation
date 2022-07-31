import os
import sys
import numpy as np
from analysis.analysis_pipeline import pipeline_replicates
from msa.tools.check_length import check_length

sysid = sys.argv[1].strip(".fa")
root = os.getcwd()
print(f"System ID: {sysid}\nRoot folder: {root}\noutput folder: {out}")

# sysid = "1a3aA"
msa = os.path.join("systems", sysid, f"{sysid}_filtered_25.fasta")
seq_l = check_length(msa)
zbool = True
dca_mode = "mf"
if zbool:
    threshold_values = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]
else:
    threshold_values = np.arange(10, 110, 10)

ppl = pipeline_replicates(sysid, seq_l, threshold_values, npairs=0, zfilter=True, plots=True, passthrough=False)
