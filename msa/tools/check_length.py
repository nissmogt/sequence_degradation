
def check_length(filein):
    from Bio import AlignIO
    align = AlignIO.read(open(filein), "fasta")
    return align.get_alignment_length()