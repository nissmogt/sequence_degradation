import numpy as np


def cm_make(score_matrix, score):
    """
    Makes a contact map from a DCA score matrix file. Also maps dca indices to pdb.
    Rank by scores.
    :param score_matrix: Frobenius norm matrix
    :param score:
    :return: Three-column dataframe composed of pair i, j, and fn score
    """
    import pandas as pd
    L = score_matrix.shape[0]
    dca_scores = []
    for i in range(L - 1):
        for j in range(i + 1, L):
            dca_scores.append([int(i+1), int(j+1), score_matrix[i, j]])    # added 1 to indices since they start from 0
    if score == 'fn':
        df_dca = pd.DataFrame(np.array(dca_scores), columns=["i", "j", score])
    if score == 'di':
        df_dca = pd.DataFrame(np.array(dca_scores), columns=["i", "j", score])
    if score == 'fn_apc':
        df_dca = pd.DataFrame(np.array(dca_scores), columns=["i", "j", score])
    if score == 'di_apc':
        df_dca = pd.DataFrame(np.array(dca_scores), columns=["i", "j", score])
    return df_dca


def calculate_fni_score(msa_name, j_matrix, j_matrix_null):
    # Function that calculates Frobenius Norm score of coupling matrices.
    q = j_matrix.shape[0]
    L = j_matrix_null.shape[1]
    fn_1 = np.zeros((L, L))
    fn_2 = np.zeros((L, L))
    fni_scores = np.zeros((L, L))
    print("Size %d" % L)

    for i in range(L - 1):
        for j in range(i + 1, L):
            fn_1[i, j] = np.linalg.norm((j_matrix[i, j]), 'fro')
            fn_2[i, j] = np.linalg.norm((j_matrix_null[i, j]), 'fro')
            # Take the difference between paired and scrambled FN
            fni_scores[i, j] = (fn_1[i, j] - fn_2[i, j])
            fni_scores[j, i] = fni_scores[i, j]

    L = fni_scores.shape[0]
    if L != fni_scores.shape[1]:
        raise ValueError("Input matrix is not symmetric: {}".format(fni_scores.shape))
    col_mean = np.mean(fni_scores, axis=0) * L / (L - 1)
    row_mean = np.mean(fni_scores, axis=1) * L / (L - 1)
    matrix_mean = np.mean(fni_scores) * L / (L - 1)

    # APC correction
    # corrnorms = fni_scores - np.outer(col_mean, row_mean) / matrix_mean
    apc = np.dot(col_mean.reshape(L, 1), col_mean.reshape(1, L)) / matrix_mean
    corrnorms = fni_scores - apc
    corrnorms[np.diag_indices(L)] = 0
    df_scores_apc = cm_make(msa_name, corrnorms)
    df_scores = cm_make(msa_name, fni_scores, apc=False)
    return df_scores, df_scores_apc


def read_matlab_matrix(filein):
    # Function that reads in the couplings and fields Matlab matrix file.
    import h5py
    mat = {}
    f = h5py.File(filein, 'r')
    for k, v in f.items():
        mat[k] = np.array(v)
    h = mat['h']
    J = mat['J']
    q = mat['h'].shape[0]
    N = mat['h'].shape[1]
    couplings = J.T
    fields = h.T
    f.close()
    return fields, couplings


def average_jmatrix(msa_name, nReplicates):
    # Function that averages over coupling matrices.
    # initialize coupling sum to zero
    sum_couplings = []
    resultDir = "scrambled_results\\"
    print("(average_jmatrix)")
    for i in range(nReplicates):
        dca_matrix = "{}matrix_ising_{}_rep{}_scrambled.fas.mat".format(resultDir, msa_name, i)
        fields, couplings = read_matlab_matrix(dca_matrix)
        sum_couplings.append(couplings)
    average_couplings = sum(sum_couplings) / nReplicates
    outfile = "scrambled_results\\{}_avg_rep{}_matrix.npy".format(msa_name, nReplicates)
    np.save(outfile, average_couplings)
    return average_couplings

pdb = "TAZ1_CITED2_hmmer"
# pdb = "TAZ1_HIF1A_hmmer"
mf = "matrix_ising_{}.fas.mat".format(pdb)
h, J = read_matlab_matrix(mf)
a = np.zeros((150,150))
ast = np.zeros((150,150))
for ri in range(149):
    for rj in range(150):
        l = []
        ll = []
        if ri != rj:
            jij_values = J[ri][rj]
            a[ri][rj] = np.mean(jij_values)
            ast[ri][rj] = np.std(jij_values)

            # for i in range(21):
            #     x = J[ri][rj][i]
                # for j in range(21):
                    # y = J[ri][rj][i][j]
                    # if min(x) < y:
                    #     l.append([y, (ri, rj)])
                # ll.append(min(x))
                    # else:
                    #     l.append([min(x), (ri, rj)])
                    #     ll.append(min(x))
            # print(ll)
            # aa = np.array(ll)

            # msa = "1GL2_A_1GL2_D"
            # vanilla_dca_matrix = "results\\matrix_files\\matrix_ising_{}.fas.mat".format(msa)
            # avg_J = average_jmatrix(msa, 5)
            # h_paired, J_paired = read_matlab_matrix(vanilla_dca_matrix)
            # df_fni, df_fni_apc = calculate_fni_score(msa, J_paired, avg_J)
