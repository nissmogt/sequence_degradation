#!/home/kmm5/.conda/envs/bio/bin/python3.8
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=fail
#SBATCH --mail-user=kmm5@rice.edu
#SBATCH --profile=ALL
#SBATCH --time=00:01:00
#SBATCH --array=0

import os
import subprocess
# Activate conda environment
#subprocess.call(["source $(conda info --base)/etc/profile.d/conda.sh \
#                  conda activate bio"], shell=True)

msa_files = os.listdir("aln")
nfiles = len(msa_files)
hash_table = dict(zip(range(nfiles), msa_files))
x = int(os.getenv("SLURM_ARRAY_TASK_ID"))
msa = (hash_table.get(x)).strip(".fa")
print(f"{msa}")
cmd=f'sbatch --job-name={msa} \
             --time=24:00:00 \
             -o logs/{msa}_{x}.out \
	     -e logs/{msa}_{x}.err \
      	     --account=ctbp-onuchic \
	     --partition=ctbp-onuchic \
	     --nodes=1 \
	     --ntasks=1 \
      	     --cpus-per-task=16 \
  	     --mem-per-cpu=8G \
	     --mail-type=fail \
	     --mail-user=kmm5@rice.edu \
      	     --profile=ALL \
      	     --wrap="python run_pipeline.py {msa}"'
subprocess.call([cmd], shell=True)
