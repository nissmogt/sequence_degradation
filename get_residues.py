def three2one(sequence):
    """ Lookup table - translate a protein sequence from 3 to 1 letter code
    """

    code = {"GLY": "G", "ALA": "A", "LEU": "L", "ILE": "I",
            "ARG": "R", "LYS": "K", "MET": "M", "CYS": "C",
            "TYR": "Y", "THR": "T", "PRO": "P", "SER": "S",
            "TRP": "W", "ASP": "D", "GLU": "E", "ASN": "N",
            "GLN": "Q", "PHE": "F", "HIS": "H", "VAL": "V",
            "M3L": "K", "MSE": "M", "CAS": "C"}

    newprot = ""
    for a in sequence:
        newprot += code.get(a, "?")

    return newprot


def get_residues(chains, seq=False):
    """
    Build a simple list of residues from a single chain of a PDB file.
    :param chains: Bio PDB Chain Object
    :param seq: Boolean (Default: False) - Outputs sequence if True.
    :return: A list of Bio.PDB.Residue objects.
    """
    import Bio.PDB
    import sys
    import os
    # from parse_pdb import parse_pdb
    # chains = parse_pdb(pdb)
    chain_lengths = []
    residues = []
    sequence = []
    for ch in chains:
        # make sure res are standard AA
        num_residues = 0
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            # if Bio.PDB.is_aa(res, standard=True):
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            num_residues += 1
            if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS' or res_id == 'HMS') and is_regular_res:
                residues.append(res)
                sequence.append(res.get_resname())
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))
        chain_lengths.append(num_residues)

    if seq:
        sequence = three2one(sequence)
        seq_a = sequence[:chain_lengths[0]]
        seq_b = sequence[chain_lengths[0]:]
        return seq_a, seq_b
    else:
        return residues, chain_lengths


def get_residues_web(pdb, seq=False, chain_ids=None):
    """
    Build a simple list of residues from a single chain of a PDB file.
    :param pdb: String - PDB filename
    :param seq: Boolean (Default: False) - Outputs sequence if True.
    :param chain_ids:
    :return: A list of Bio.PDB.Residue objects.
    """
    import Bio.PDB
    import sys
    import os

    if chain_ids is None:
        chain_ids = ["A", "B"]

    # chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]
    chain_lengths = []
    parser = Bio.PDB.PDBParser()

    struct = parser.get_structure(pdb, pdb)
    model = struct[0]
    # if len(self.chain_ids) == 0:
    # get residues from every chain.
    #    chains = model.get_list()
    # else:
    chains = [model[ch_id] for ch_id in chain_ids]

    residues = []
    sequence = []
    for ch in chains:
        # make sure res are standard AA
        num_residues = 0
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            # if Bio.PDB.is_aa(res, standard=True):
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
                residues.append(res)
                sequence.append(res.get_resname())
                num_residues += 1
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))
        chain_lengths.append(num_residues)

    if seq:
        sequence = three2one(sequence)
        seq_a = sequence[:chain_lengths[0]]
        seq_b = sequence[chain_lengths[0]:]
        return seq_a, seq_b
    else:
        return residues, chain_lengths
