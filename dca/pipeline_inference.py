import os
import matlab.engine

usr = os.path.expanduser("~")
#dir_dca = os.path.join(usr, "MATLAB-Drive", "Research", "Scripts", "mf_plm_reweight")


def inference(_pfamid, dir_dca, matlab_input, model_length, pseudocount=0.2):
    output_dir = os.path.join(os.path.dirname(matlab_input), "mf", f"pc{pseudocount:.1f}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_dca_output = os.path.join(output_dir, f"DI_{_pfamid}_n{model_length}.txt")
    file_matrix_output = os.path.join(output_dir, f"matrix_{_pfamid}_n{model_length}.mat")
    if not os.path.exists(file_dca_output):
        eng = matlab.engine.start_matlab()
        eng.addpath(dir_dca, nargout=0)
        eng.addpath(os.path.join(dir_dca, "functions"))
        eng.addpath(output_dir, nargout=0)
        neffective = eng.dca_h_J_Full_v5(0.2, matlab_input, file_dca_output, file_matrix_output, nargout=1)
        eng.quit()
        return neffective
    else:
        return 0
