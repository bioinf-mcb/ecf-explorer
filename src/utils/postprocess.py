import pandas as pd

def load_clustering_results(workdir_path):
    """Load the results of clustering of hhblits hits to proteins. Input file
    should be stored in csv format.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.

    Returns:
    --------
    pd.DataFrame
        Pandas dataframe with loaded data.
    """
    results_filepath  = workdir_path \
    + '/ecf-explorer/output/reprseq-all-clust.txt'
    clustering_results = pd.read_csv(results_filepath, sep=',')

    return clustering_results

def select_reprseqs_to_display(clustering_results, reprseqs_to_display):
    """From the dataframe of hhblits hits to proteins, select only those rows
    that contains data from representative sequences with the numbers given in
    reprseqs_to_display list.

    Parameters:
    -----------
    clustering_results : pd.DataFrame
        Pandas dataframe with results of all vs all hhblits search of
        representative proteins that was already processed with clustering.
    reprseqs_to_display : list
        List of numerical identifiers of representative proteins. The input
        dataframe will be filtered so it will contain data regarding only those
        proteins.

    Returns:
    --------
    pd.DataFrame
        Pandas dataframe with filtered data.
    """

    # convert list of numerical identifiers to full string identifiers
    rseq_id_width = int(clustering_results.iloc[0].str.count(r'\d')['qname'])
    reprseq_ids = ['reprseq{:0>{idw}}'.format(rseq, idw=rseq_id_width) \
    for rseq in reprseqs_to_display]

    # select data from input dataframe
    select_reprseq_data = \
    clustering_results[clustering_results['qname'].isin(reprseq_ids)]

    return select_reprseq_data
