import pandas as pd
import glob
from subprocess import call

def split_dataset_df_into_single_proteins(workdir_path):
    """Split dataframe with filtered results of HHsuite scan into separate
    files, each file one protein.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.

    Returns:
    --------
    """
    # clean save directory
    dir_path_save = workdir_path + '/ecf-explorer/intermediate/protein_dfs/'
    if len(glob.glob(dir_path_save + '*.*')) > 0:
        call('rm  ' + dir_path_save + '*.*', shell=True)
    # load data
    dataset_df = pd.read_csv(workdir_path \
    + '/ecf-explorer/intermediate/table-hhr_filtered.txt')
    # groupby reprseq and save
    for prot_id, prot_df in dataset_df.groupby('qname'):
        prot_df.to_csv(dir_path_save + prot_id + '.txt', index=False)

def create_bash_script_to_run_clustering(workdir_path, ecf_env_path,
                                         lib_ecf_path, dbscan_params):
    """Create bash script to run clustering on all proteins in the dataset in
    a parallel manner.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.
    dbscan_params : dict
        Dictionary of parameters to be used in DBSCAN. Two parameters must be
        set (see scikit-learn docs for descriptions):
        'eps' <float>
        'min_samples' <int>
    ecf_env_path : string
        Path to virtualenv where ecf-explorer code is stored. It needs to be
        communicated to the script to activate virtualenv in bash script.
    lib_ecf_path : string
        Path to ecf-explorer code. It needs to be communicated to the Python
        script so the clustering functions will be loaded.

    Returns:
    --------
    string
        Absolute path to the created bash script.
    """
    script_filepath = workdir_path + '/ecf-explorer/scripts/' \
    + 'run_clustering.sh'

    cmd = '#!/bin/bash\n\n'
    cmd += 'FILE=$(basename "${1}")\nFILE=${FILE%.*}\n'
    cmd += 'source ' + ecf_env_path + '/bin/activate\n'
    cmd += 'python3 ' + ecf_env_path \
    + '/ecf-explorer/src/ecfs/run_clustering_single_protein.py ${1} ' \
    + str(dbscan_params['eps']) + ' ' + str(dbscan_params['min_samples']) \
    + ' ' + lib_ecf_path + '\n'

    with open(script_filepath, 'w') as file_sh:
        file_sh.write(cmd)

    return script_filepath

def run_clustering_with_bash(workdir_path, n_cores):
    """Run previously prepared bash script that will run clustering on all
    proteins in dataset in a parallel manner.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.
    n_cores : int
        Maximum number of cores to be used in parallel run.

    Returns:
    --------
    """
    script_filepath = workdir_path + '/ecf-explorer/scripts/' \
    + 'run_clustering.sh'
    protein_df_dirpath = workdir_path + '/ecf-explorer/intermediate/protein_dfs'

    # set execute mode
    call('chmod a+x ' + script_filepath, shell=True)

    # run script
    run_cmd = 'nohup find {} -name "reprseq*.txt" | xargs -P {} -n 1 {} \
    &'.format(protein_df_dirpath, n_cores, script_filepath)
    call(run_cmd, shell=True)

def concatenate_clustering_results(workdir_path):
    """Check if all proteins were processed with clustering and if yes then
    colect all the resulting dataframes into summary files. Separate file will
    be created for clustering results files and cluster statistics files. If
    not all expected output files are present function quits.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.

    Returns:
    --------
    bool
        Status of the execution. Returns False if still waiting for processing
        to complete. Returns True if process complete and files were successfuly
        concatenated.
    """
    # check if all expected output files are present
    protein_df_dirpath = workdir_path + '/ecf-explorer/intermediate/protein_dfs'
    nb_expected_files  = len(glob.glob(protein_df_dirpath + '/*.txt'))
    nb_complete_files  = len(glob.glob(protein_df_dirpath + '/*.out'))

    if 2 * nb_expected_files == nb_complete_files:
        print('All files processed. Concatenating...')
        proteins_clust = []
        for prot_df_clust in glob.glob(protein_df_dirpath + '/*.clust.out'):
            proteins_clust.append(pd.read_csv(prot_df_clust))
        proteins_clust_all = pd.concat(proteins_clust)
        proteins_clust_all.to_csv(workdir_path + '/ecf-explorer/output/' \
                                  + 'reprseq-all-clust.txt', index=False)
        proteins_stat = []
        for prot_df_stat in glob.glob(protein_df_dirpath + '/*.stat.out'):
            proteins_stat.append(pd.read_csv(prot_df_stat))
        proteins_stat_all = pd.concat(proteins_stat)
        proteins_stat_all.to_csv(workdir_path + '/ecf-explorer/output/' \
                                  + 'reprseq-all-stat.txt', index=False)
        print('Done.')
        return True
    else:
        print('Not all files processed. Check again later.')
        return False


def clean_clustering_partial_data(workdir_path):
    """Clean directory with single protein dataframes. Saves space.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.

    Returns:
    --------
    """
    dir_path_single = workdir_path + '/ecf-explorer/intermediate/protein_dfs/'
    if len(glob.glob(dir_path_single + '*.*')) > 0:
        call('rm  ' + dir_path_single + '*.*', shell=True)
    print('Cleaned single protein dataframes.')
