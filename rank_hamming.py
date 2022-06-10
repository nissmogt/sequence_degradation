def rank_hamming(df, score, threshold=5):
    """
    Rank DCA predictions by designated score and remove pairs with sequence distance < threshold
    :param df: Dataframe with four columns i,j, fn_apc, and fn
    :param score: String, which column to sort by
    :param threshold: int, Distance threshold
    :return: Dataframe ranked by score and filtered
    """
    df_sorted = df.sort_values(ascending=False, by=[score])
    df_hamming = df_sorted[abs(df_sorted['i'] - df_sorted['j']) >= threshold].reset_index(drop=True)
    return df_hamming
