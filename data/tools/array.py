import numpy as np


def save_array(array, outname):
    np.save(outname, array)


def calculate_ratio_below_threshold(array, threshold):
    return np.count_nonzero(array <= threshold) / len(array)