import os
import analysis.zscore
import pandas as pd


class DirectCoupling:
    def __init__(self):
        self._sysid = ""
        self._pdbid = ""
        self._score = ""
        self._template_sequence = ""

    def set_sysid(self, sysid):
        self._sysid = sysid

    def set_pdbid(self, pdbid):
        self._pdbid = pdbid

    def set_score(self, score):
        self._score = score

    def set_template_sequence(self, template):
        self._template_sequence = template

    def load_to_df(self, dca_in):
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

    def add_pdb_distances(self, dca_dataframe, pdb_dataframe):
        """
        Uses pandas merge to get the intersection between dca and pdb pairs
        """
        return dca_dataframe.merge(pdb_dataframe, how='inner', on=['i', 'j'])

    def zscore(self, dca_dataframe, reference):
        return analysis.zscore.zscore_calculation(dca_dataframe, self._score, reference)

    def index_shift(self, dca_dataframe, cols, shift):
        """
        cols: tuple of strings
        shift: integer

        """
        c1, c2 = cols
        dca_dataframe[c1] = dca_dataframe[c1] + shift
        dca_dataframe[c2] = dca_dataframe[c2] + shift
        return dca_dataframe
