#!/home/kmm5/.conda/envs/bio/bin/python3.8
import os
import sys
import data.system_object

sysid = sys.argv[1].strip(".fa")
aln_dir = sys.argv[2]
root = os.getcwd()
out = "/scratch/kmm5/"
matlab_dir = os.path.join(root, "mf_plm_reweight")
print(f"System ID: {sysid}\nRoot folder: {root}\noutput folder: {out}")

s = data.system_object.System()
s.set_sysid(sysid)
s.set_align_dir(os.path.join(root, aln_dir))
s.make_new_dirs(root, out)
s.check_dir()
s.filter(run=True)
list_len = s.replicates(run=True)
list_neff = s.run_inference(_nreplicates=70, _dir_dca=matlab_dir, _len_list=list_len, passthrough=False)
