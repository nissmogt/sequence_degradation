import os
import numpy as np
import msa.tools.gap_filter
import msa.msa_functions
import dca.pipeline_inference


class System:
    def __init__(self):
        self._sysid = ""
        self._dir_aln = ""

        self._dir_sys = ""
        self._dir_results = ""
        self._dir_replicates = ""
        self._dir_avg_results = ""
        self._raw_msa = ""
        self._filtered_msa = ""
        self._nseq = 0

    def set_sysid(self, sysid):
        self._sysid = sysid

    def set_align_dir(self, directory):
        self._dir_aln = directory

    def make_new_dirs(self, _root, _out, replicates=True):
        # DIRECTORY STRUCTURE INITIALIZATION
        self._dir_sys = os.path.join(_out, "systems", self._sysid)
        self._dir_results = os.path.join(self._dir_sys, "results")
        self._dir_avg_results = os.path.join(self._dir_results, "average_ppv")
        if replicates:
            self._dir_replicates = os.path.join(self._dir_sys, "replicates")
            list_dir = [self._dir_sys, self._dir_results, self._dir_avg_results, self._dir_replicates]
        else:
            list_dir = [self._dir_sys, self._dir_results, self._dir_avg_results]

        # MAKE DIRECTORY IF DOESNT EXIST
        for entry in list_dir:
            if not os.path.exists(entry):
                os.makedirs(entry)

    def check_dir(self):
        print(f"SYS: {self._dir_sys}, RESULTS: {self._dir_results}, "
              f"AVG_RES: {self._dir_avg_results}, REPLICATES: {self._dir_replicates}")

    def filter(self, pct_gaps=0.2, run=True):
        """
        Applies 25% gap filter to msa.
        """
        self._raw_msa = os.path.join(self._dir_aln, f"{self._sysid}.fa")

        if os.path.exists(self._raw_msa):
            self._filtered_msa = os.path.join(self._dir_sys, f"{self._sysid}_full_filtered_25.fasta")
        else:
            print(f"{self._raw_msa} does not exist.")
        if run:
            self._filtered_msa, self._nseq, percent_gaps = msa.tools.gap_filter.gap_filter(self._raw_msa,
                                                                                           pct_gaps, self._dir_sys)

    def replicates(self, run=True):
        """
        Read list of lengths generated by gen_replicates
        """
        nseq_list = os.path.join(self._dir_replicates, "length_list.txt")
        if not run:
            if os.path.exists(nseq_list):
                with open(nseq_list, "r") as fp:
                    list_of_len = [int(line.rstrip()) for line in fp]
                    assert list_of_len[-1] > 0
        if run:
            print("Running replicate code")
            # 3. Replicate generation 100 replicates
            list_of_len = msa.msa_functions.generate_replicates(self._filtered_msa, 100, self._dir_replicates)
        return list_of_len

    def run_inference(self, _nreplicates, _dir_dca, _len_list=None, passthrough=False):
        """
        Infer DCA couplings for every downsampled ensemble of replicates
        """
        if _nreplicates == 1:
            # For one DCA run
            output = os.path.join(self._dir_sys, "neff_array.npy")
            msa_input = self._filtered_msa
            model_length = self._nseq
            n_effective = dca.pipeline_inference.inference(self._sysid, _dir_dca, msa_input, model_length)
            np.save(output, n_effective)
            return n_effective

        else:
            output = os.path.join(self._dir_replicates, "neff_array.npy")
            if passthrough:
                # Load neffective array. FUTURE: Begin from last replicate
                return np.load(output)
            else:
                nmodels = len(_len_list)
                n_effective_array = np.zeros((_nreplicates, nmodels))  # initialize the array by number of replicates
                for rep in range(_nreplicates):
                    msa_rep_dir = os.path.join(self._dir_replicates, f"sub{rep}")
                    for model in range(nmodels):
                        model_length = _len_list[model]
                        if model_length > 0:
                            msa_input = os.path.join(msa_rep_dir, f"{self._sysid}_n{model_length}_sub{rep}.fasta")
                            print(f"PFAM: {self._sysid} REP: {rep} N{model}: {model_length}")
                            print(f"{msa_input}")
                            n_effective_array[rep][model] = dca.pipeline_inference.inference(self._sysid, _dir_dca,
                                                                                             msa_input, model_length)
                    # Save Neff to file after every replicate
                    np.save(output, n_effective_array)
                return n_effective_array
