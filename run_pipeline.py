import os
import data.system_object

pfamid = "1a6mA"
root = "/Users/euler/PycharmProjects/sequence_degradation/"

sys = data.system_object.System()
sys._sysid = pfamid
sys._dir_aln = os.path.join(root, "aln")
sys._dir_pdb = os.path.join(root, "pdb")
sys.make_new_dirs(root)
sys.filter()
sys.replicates()
