def parse_pdb(pdb_in):
    import Bio.PDB
    import os
    if os.path.dirname(pdb_in):
        dir_pdb = os.path.dirname(pdb_in)
        pdb_file = os.path.basename(pdb_in)
    else:
        dir_pdb = ""
        pdb_file = pdb_in

    pdb_id, ext = os.path.splitext(pdb_file)
    chain_lengths = []
    parser = Bio.PDB.PDBParser()
    struct = parser.get_structure(pdb_id, pdb_in)
    model = struct[0]
    chains = [ch_id for ch_id in model.get_list()]
    return chains