def calculate_ppv(dataframe_in, threshold):
    import numpy as np

    dca_dataframe = dataframe_in.reset_index(drop=True)
    tp_count = np.zeros(len(dca_dataframe))
    fp_count = np.zeros_like(tp_count)
    positive_predictive_value = np.zeros_like(tp_count)
    tcount = 0
    fcount = 0

    for j in range(0, len(dca_dataframe)):
        if dca_dataframe["d"][j] <= threshold:
            tcount += 1
            tp_count[j] = tcount
            fp_count[j] = fcount
        else:
            fcount += 1
            tp_count[j] = tcount
            fp_count[j] = fcount
        positive_predictive_value[j] = tp_count[j] / (fp_count[j] + tp_count[j])

    return positive_predictive_value, tp_count + fp_count
