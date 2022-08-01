#!/home/kmm5/.conda/envs/bio/bin/python3.8
#SBATCH --job-name=an_1a3aA
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
#SBATCH -o logs/a_1a3aA.out
#SBATCH -e logs/a_1a3aA.err

# conda activate bio
python run_analysis.py
