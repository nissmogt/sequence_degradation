import os
import analysis.zscore
import pandas as pd


class DirectCoupling:
    def __init__(self):
        self._sysid = ""
        self._pdbid = ""
        self._score = ""
        self._template_sequence = ""

    def load_raw(self, dca_in):
        df_header = ["i", "j", "di", "zscore", "diapc", "mi", "d", "si", "sj", "chain_1", "chain_2", "resnames",
                     "atom_id"]
        if os.path.exists(dca_in):
            dir_dca_in = os.path.dirname(dca_in)
            return pd.read_csv(dca_in, delimiter=',', header=0)
        else:
            print(f"File {dca_in} doesn't exist.")

    def savetocsv(self, dataframe, _filename):
        import pandas as pd
        df_header = ["i", "j", "di", "zscore", "diapc", "mi", "d", "si", "sj", "chain_1", "chain_2", "resnames",
                     "atom_id"]
        dataframe.to_csv(_filename, sep='\t', index=False, header=df_header, float_format='%.4f')

    def rank_hamming(self, dca_dataframe, distance=5):
        """
        Rank DCA predictions by designated score and remove pairs with sequence distance < threshold
        :param dca_dataframe: Dataframe
        :param distance: int, Distance threshold
        :return: Dataframe ranked by score and filtered
        """
        df_sorted = dca_dataframe.sort_values(ascending=False, by=[self._score])
        df_ranked = df_sorted[abs(df_sorted['i'] - df_sorted['j']) >= distance].reset_index(drop=True)
        return df_ranked

    def zscore(self, dca_dataframe, reference):
        return analysis.zscore.zscore_calculation(dca_dataframe, self._score, reference)
