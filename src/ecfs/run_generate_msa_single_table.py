import sys
import pandas as pd

if __name__ == '__main__':
    workdir_path = sys.argv[1]
    i            = sys.argv[2]
    ident_cut    = float(sys.argv[3])
    qcov_cut     = float(sys.argv[4])
    eval_cut     = float(sys.argv[5])

    sys.path.append(sys.argv[6])
    from utils.postprocess import load_clustering_results
    from ecfs.msa import save_qecf_as_msa

    parts_dir = workdir_path + '/ecf-explorer/intermediate/qecf_dfs/'
    clustering_stats = pd.read_csv(parts_dir + 'clustering_stats-' + str(i) \
                                   + '.csv')
    clustering_results = load_clustering_results(workdir_path, mode='clusters')

    msa_params = {'ident_cut':ident_cut, 'qcov_cut':qcov_cut,
                  'eval_cut':eval_cut}

    save_qecf_as_msa(workdir_path, clustering_stats, clustering_results,
                     msa_params)
