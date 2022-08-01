import os
import sys
import data.system_object

sysid = sys.argv[1].strip(".fa")
root = os.getcwd()
out = "/scratch/kmm5/"
matlab_dir = os.path.join(root, "mf_plm_reweight")
print(f"System ID: {sysid}\nRoot folder: {root}\noutput folder: {out}")

s = data.system_object.System()
s._sysid = sysid
s._dir_aln = os.path.join(root, "aln")
s._dir_pdb = os.path.join(root, "pdb")
s.make_new_dirs(root, out)
s.filter()
list_len = s.replicates()
list_neff = s.run_inference(list_len, 100, matlab_dir, passthrough=False)