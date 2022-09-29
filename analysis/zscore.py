import numpy as np


def add_systems(scores, reference):
    return np.concatenate((reference, scores))


def zscore_calculation(dca_dataframe, dca_score, reference):
    # calc zscore
    # reference = np.load(f"dimers_{dca_score.upper()}_scores.npy")
    scores = dca_dataframe[dca_score].to_numpy()
    mean = np.mean(reference)
    std = np.std(reference)
    z = np.zeros_like(scores)
    for k in range(len(scores)):
        z[k] = (scores[k] - mean) / std
    dca_dataframe.insert(3, "zscore", z)
    return dca_dataframe


def calculate_average_zscore(array):
    return np.mean(array, axis=0), np.std(array, axis=0)


def make_reference(_list, score_type, outf):
    ref = np.concatenate(_list)
    np.save(outf, ref)
    return ref

