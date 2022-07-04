import numpy as np


def zscore_calculation(dca_dataframe, dca_score, reference):
    # calc zscore
    # reference = np.load(f"dimers_{dca_score.upper()}_scores.npy")
    if not reference:
        reference = np.load(f"monomer_{dca_score.upper()}_scores.npy")
    else:
        scores = dca_dataframe[dca_score].to_numpy()
        new_reference = np.concatenate((reference, scores))
        mean = np.mean(new_reference)
        std = np.std(new_reference)
        z = np.zeros_like(scores)
        for k in range(len(scores)):
            z[k] = (scores[k] - mean) / std
        dca_dataframe.insert(3, "zscore", z)
        return dca_dataframe
