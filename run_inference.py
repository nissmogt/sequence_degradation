#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import sys
import data.system_object

# Single DCA run
sysid = sys.argv[1].strip(".fa")
aln_dir = sys.argv[2]
# sysid = "1cc8A"
# aln_dir = "aln"
root = os.getcwd()
out = "/scratch/kmm5/single/"
# out = "tests/single/"
if not os.path.exists(out):
    os.makedirs(out)
matlab_dir = os.path.join(root, "mf_plm_reweight")
print(f"System ID: {sysid}\nRoot folder: {root}\noutput folder: {out}")

s = data.system_object.System()
s.set_sysid(sysid)
s.set_align_dir(os.path.join(root, aln_dir))
s.make_new_dirs(root, out, replicates=False)
s.check_dir()
# for a range of pct_gaps [0.05, 0.1, 0.15, 0.2]
gap_list = [0.05, 0.1, 0.15, 0.2]
s.filter(run=True, pct_gaps=gap_list[1])
neff = s.run_inference(_len_list=None, _nreplicates=1, _dir_dca=matlab_dir, passthrough=False)
print(neff)
