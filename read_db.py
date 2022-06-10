import numpy as np
df_file = "C:\\Users\\kmehr\\PycharmProjects\\idp_trimer\\media-3_added.xlsx"


def get_lengths(msa_name):
    import pandas as pd
    df_db = pd.read_excel(df_file)

    # used in query() pandas function
    pdbid = msa_name
    chain_1 = 'A'
    chain_2 = 'B'

    try:
        x = df_db.query("PDB == @pdbid and Chain_1 == @chain_1 and Chain_2 == @chain_2")
        chain_1_length = (x.Sequence_Length_1.values[0]).astype('int64')
        chain_2_length = (x.Sequence_Length_2.values[0]).astype('int64')
        return [chain_1_length, chain_2_length]
    except ValueError:
        print("{} not found in query database.".format(msa_name))
        return [1, 1]
