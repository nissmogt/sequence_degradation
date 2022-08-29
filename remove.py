#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import subprocess
# Activate conda environment
#subprocess.call(["source $(conda info --base)/etc/profile.d/conda.sh \
#                  conda activate bio"], shell=True)

for i in range(102):
	msa_files = os.listdir("aln")
	nfiles = len(msa_files)
	hash_table = dict(zip(range(nfiles), msa_files))
	x = i
	msa = (hash_table.get(x)).strip(".fa")
	print(f"{msa}")
	for repid in range(100):
		if os.path.exists(f"/scratch/kmm5/systems/{msa}/replicates/sub{repid}/mf/pc0.2/images"):
			print(f"/scratch/kmm5/systems/{msa}/replicates/sub{repid}/mf/pc0.2/images")
			cmd=f'rm -r /scratch/kmm5/systems/{msa}/replicates/sub{repid}/mf/pc0.2/images'
			subprocess.call([cmd], shell=True)
