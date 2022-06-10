temp_mdlnum.txt]
overall how many models are analyzed for all-atom-based ensemble

temp_atmnum.txt]
Total number of atoms being analyzed

temp_fitted_model.pdb]
all-atom-based best-fitted ensemble through an iterative process 

temp_avg.pdb]
the average fine-grained structure (all atoms)

temp_avg_close_model_no.txt]
the model that is closest to temp_avg.pdb (having the smallest rmsd)

temp_avg_close_model.pdb]
coordinates of the model that is closest to temp_avg.pdb (having the smallest rmsd)

temp_rmsd.txt]
root mean square fluctuation of every atom

temp_avg_rmsd.txt]
average root mean square fluctuations per atom

temp_eg_values.txt]
eigenvalues obtained from singular value decomposition of all-atom-based covariance

temp_egvals.txt]
eigenvalues obtained from eigen value decomposition of all-atom-based covariance (the square of values in temp_eg_values.txt)

temp_eg_vectors.txt]
eigenvectors obtained from singular value decomposition of all-atom-based covariance

temp_trj.txt]
PC ccodinates of each model on different modes (all-atom-based)

temp_effect_freqs.txt]
mode frequencies at 300K based on the variance of the modes (eigenvalues at temp_egvals.txt) and harmonic approximation





temp_ca_mdlnum.txt]
overall how many models are analyzed for node-based ensemble

temp_canum.txt]
number of nodes being analyzed

temp_fitted_model_ca.pdb]
node-based best-fitted ensemble through an iterative process 

temp_ca_avg.pdb]
the average coarse-grained structure (nodes only)

temp_avg_close_model_no_ca.txt]
the model that is closest to temp_ca_avg.pdb (having the smallest rmsd)

temp_avg_close_model_ca.pdb]
coordinates of the model that is closest to temp_ca_avg.pdb (having the smallest rmsd)

temp_ca_rmsd.txt] & temp.rmsd]
root mean square fluctuation of every node

temp_ca_avg_rmsd.txt]
average root mean square fluctuations per node

temp_ca_eg_values.txt]
eigenvalues obtained from singular value decomposition of node-based covariance

temp_ca_egvals.txt]
eigenvalues obtained from eigen value decomposition of node-based covariance (the square of values in temp_ca_eg_values.txt)

temp_ca_eg_vectors.txt]
eigenvectors obtained from singular value decomposition of node-based covariance

temp_ca_trj.txt]
PC ccodinates of each model on different modes (node-based)




temp_mdl_ori.txt]
how many models are in the originally input file

temp.chain]
which chain being selected

temp.cont]
number of C-alphas in contact within a cutoff distance of 10 Angstrom (based on the first model in the temp_fitted_model_ca.pdb)

temp.nodes]
the taken representative atoms from amino acid (CA) and nucleotides (Phosphorus/C4*(in the sugar)/C2 (in the base))

temp.sloweigenvectors]
same as 'temp_ca_eg_vectors.txt'

temp.slowmodes]
magnitude square of each node on every mode

temp.slowav]
weighted magnitude square of each node on every mode from the first two modes in 'temp.slowmodes'
