
def vectorize_pdb_contacts(pdb_df, dimer_length):
    """
    Converts PDB pairs into binary matrix and then flattens into 1-D array
    :param pdb_df: Dataframe object - PDB pairs
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    import pandas
    import numpy as np
    pdb_array = pdb_df.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Initialize contact matrices
    pdb_matrix = np.zeros((dimer_length, dimer_length))

    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in pdb_array:
        pdb_matrix[int(i) - 1, int(j) - 1] = 1

    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    pdb_flat_matrix = pdb_matrix.ravel()
    return pdb_flat_matrix
