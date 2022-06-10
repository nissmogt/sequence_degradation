using ParalogMatching

# infile2 = "cited2_hmmer_1r8u_filtered_25.fasta"
infile2 = "hif1a_1l8u_filtered_25.fasta"
infile1 = "taz1_7lvs_filtered_25.fasta"

outfile = "TAZ1_HIF1A_hmmer.fas"
X1 = read_fasta_alignment(infile1, header_regex=r"^(?<id>[^_]+)_(?<species>[^/]+)")
X2 = read_fasta_alignment(infile2, header_regex=r"^(?<id>[^_]+)_(?<species>[^/]+)")


cutoff = 500
batch = 1
strategy = "covariation"
pseudo_count = 0.5
X12 = prepare_alignments(X1, X2, cutoff=cutoff)
match = run_matching(X12, batch=batch, strategy=strategy, pseudo_count=pseudo_count)

write_fasta_match(X12, match, outfile)

println("done")
