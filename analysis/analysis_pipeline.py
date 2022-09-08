import os
import numpy as np
import data.tools.pdb as dpdb
import analysis.dca_object as dca
import analysis.plots
from analysis.validation import calculate_ppv


def process_dca(root, _sysid, _df_pdb, _nseqs, _neff, _rep, zcalc=False, shift=0):
    """
    Processes DCA output and returns a dataframe.

    """
    # File name definitions and directory creation
    raw_dca = f"DI_{_sysid}_n{_nseqs}.txt"
    dir_dca = os.path.join(root, "systems", _sysid, "replicates", f"sub{_rep}")
    dca_in = os.path.join(dir_dca, raw_dca)
    outfile = os.path.join(dir_dca, f"{_sysid}_neff{_neff}_pc0.2_all.txt")

    # Instantiate DCA object
    d = dca.DirectCoupling()
    d.set_sysid(_sysid)
    d.set_score("diapc")
    neff_array = np.load(os.path.join(root, "systems", _sysid, "replicates", "neff_array.npy"))
    n_effective = neff_array[0][0]
    out_dca = os.path.join(dir_dca, f"{_sysid}_neff{n_effective}_all.txt")

    df = d.load_to_df(dca_in)
    if shift != 0:
        df_shift = d.index_shift(df, cols=("i", "j"), shift=shift)
        df = df_shift
    df_rank = d.rank_hamming(df, distance=5)
    df_dca = d.add_pdb_distances(df_rank, _df_pdb)

    if zcalc:
        df_header = ["i", "j", "di", "zscore", "diapc", "mi", "d", "si", "sj", "chain_1", "chain_2", "resnames",
                     "atom_id"]
        zscore_reference = np.load(os.path.join("assets", f"monomer_DIAPC_scores.npy"))
        df_z = d.zscore(df_dca, zscore_reference)
        df_z.to_csv(outfile, sep='\t', index=False, header=df_header, float_format='%.4f')

        # d.savetocsv(df_z, out_dca)
        # analysis.plots.contact_map_single(df_z, monomer=df_pdb, n=10, x="test", distance_cutoff=8)
        # ppv = calculate_ppv(df_z, 6)
        return df_z
    else:
        return df_dca


def pipeline_replicates(_dca_dir, _sysid, _ncols, thresholds_list, npairs=0, zfilter=True, plots=False, passthrough=False):
    """
    For every replicate and model, process raw DCA results, calculate Zscores, and PPV
    :param _dca_dir
    :param _sysid
    :param _ncols:
    :param npairs: to use in plotting ranked pairs
    :param plots: Bool; Plots score distribution, top pairs, and zscore filtered pairs depending on zfilter.
    :param zfilter: Bool; If False, activates top pair program. If true, zscore filter program activates.
    :param passthrough: Bool; If True, skips this program and loads ppv array
    :param thresholds_list: List; List of Zscore cutoffs
    :return:
    """
    distance_cutoff = 8
    dir_system = os.path.join(_dca_dir, "systems", _sysid)
    dir_replicates = os.path.join(dir_system, "replicates")
    dir_results = os.path.join(dir_system, "results")
    dir_pdb = os.path.join(_dca_dir, "pdb")
    # Read list of lengths generated by gen_replicates
    nseq_list = os.path.join(dir_replicates, "length_list.txt")
    if os.path.exists(nseq_list):
        with open(nseq_list, "r") as fp:
            model_lengths = [int(line.rstrip()) for line in fp]
            assert model_lengths[-1] > 0
    neffective_array = np.load(os.path.join(dir_replicates, "neff_array.npy"))
    avg_neff = neffective_array.mean(axis=0)
    # CONTROL PPV OUTFILE NAME DEPENDING ON ZFILTER VALUE
    dca_score = "DIAPC"
    pdbid = _sysid
    dir_ppv = os.path.join(dir_results, "average_ppv")
    if not os.path.exists(dir_ppv):
        os.makedirs(dir_ppv)
    if zfilter:
        output_ppv = os.path.join(dir_ppv, f"{pdbid}_ppv_zscore_{dca_score}_100reps.npy")
        output_pair_rank = os.path.join(dir_ppv, f"{pdbid}_ppv_pairrank_z_{dca_score}_reps100.npy")
    else:
        if npairs > 0:
            output_ppv = os.path.join(dir_ppv, f"{pdbid}_ppv_top{npairs}_{dca_score}_100reps.npy")
        else:
            output_ppv = os.path.join(dir_ppv, f"{pdbid}_ppv_top{thresholds_list[-1]}_{dca_score}_100reps.npy")
            output_pair_rank = os.path.join(dir_ppv,
                                            f"{pdbid}_ppv_pairrank_top{thresholds_list[-1]}_{dca_score}_reps100.npy")

    # Calculate pdb distance matrix and output to file
    df_pdb = dpdb.pipeline_pdb(_sysid, dir_pdb)
    if passthrough:
        # IF ALREADY CALCULATED LOAD PPV FILE
        return np.load(output_ppv), np.load(output_pair_rank)

    else:
        # ELSE CALCULATE PPV
        n_replicates, n_sys = neffective_array.shape

        if npairs > 0:
            pos_pred_list = np.zeros((len(thresholds_list), n_replicates, n_sys, npairs))
            pair_rank_array = np.zeros_like(pos_pred_list)
        else:
            pos_pred_list = np.zeros((len(thresholds_list), n_replicates, n_sys))
            pair_rank_array = np.zeros_like(pos_pred_list)
        for rep_id in range(6):
            # Make directories for results and plots
            dir_dca_results = os.path.join(dir_replicates, f"sub{rep_id}")
            dir_contact_map = os.path.join(dir_dca_results, "images")
            if not os.path.exists(dir_contact_map):
                os.makedirs(dir_contact_map)

            for model_id in range(len(model_lengths)):
                n_degraded_seqs = avg_neff[model_id]
                raw_dca_output = os.path.join(dir_dca_results, f"DI_{_sysid}_n{n_degraded_seqs}.txt")

                neff = int(neffective_array[rep_id][model_id])
                # OUTPUT ZSCORE
                df_dca = process_dca(_dca_dir, _sysid, df_pdb, model_lengths[model_id], neff, rep_id, zcalc=True,
                                     shift=0)
                map_idx = 0
                # draw_publish_dca(_pfamid, df_dca, "A", dir_contact_map)
                # matrix = into_matrix(df_dca)
                assert df_dca is not None
                print(f"rep: {rep_id} model: {n_degraded_seqs} neff: {neff}  L: {_ncols} map:{map_idx}")
                if plots:
                    # plot_score_distribution(dca_score, n_degraded_seqs, df_dca, neff, ncols, dir_contact_map)
                    if zfilter:
                        print("zfilter")
                        analysis.plots.plot_dist_distribution(df_dca[:750], n_degraded_seqs, _ncols, dir_contact_map)
                        analysis.plots.plot_score_distribution("zscore", n_degraded_seqs, df_dca, neff, _ncols,
                                                               dir_contact_map)

                for k in range(len(thresholds_list)):
                    threshold_value = thresholds_list[k]
                    if zfilter:
                        df_filtered = df_dca[df_dca["zscore"] >= threshold_value]
                        if plots:
                            print(f"{threshold_value}")
                            # zscore_contact_map(df_filtered, pdbid, neff, ncols + map_idx, df_pdb, distance_cutoff,
                            #                    dir_contact_map, threshold_value, False)
                            analysis.plots.plot_top_zscore(df_dca, 10, n_degraded_seqs, _ncols, dir_contact_map)
                            analysis.plots.figure_1(df_dca, df_filtered, df_pdb, pdbid, distance_cutoff,
                                                n_degraded_seqs, threshold_value, _ncols, map_idx, dir_contact_map)
                    else:
                        df_filtered = df_dca[:threshold_value]

                    if len(df_filtered) > 0:
                        ppv, pair_rank = calculate_ppv(df_filtered, distance_cutoff)
                        if plots:
                            analysis.plots.plot_ppv(zfilter, ppv, neff, _sysid, _ncols,
                                                    threshold_value, dir_contact_map)
                            print("plotting PPV")
                        if npairs > 0:  # to use in plotting ranked pairs
                            pos_pred_list[k][rep_id][model_id] = ppv
                            pair_rank_array[k][rep_id][model_id] = pair_rank
                        else:
                            pos_pred_list[k][rep_id][model_id] = ppv[-1]
                            pair_rank_array[k][rep_id][model_id] = pair_rank
                    else:
                        pos_pred_list[k][rep_id][model_id] = 0
        np.save(output_ppv, pos_pred_list)
        np.save(output_pair_rank, pair_rank_array)
        return pos_pred_list, pair_rank_array
