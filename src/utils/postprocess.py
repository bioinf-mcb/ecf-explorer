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

def load_features_data(workdir_path, features_to_load):
    """Load the results of annotating repsreq with various tools. Input files
    should be stored in csv format. Each file should represent single feature.
    The csv file should have following columns:
    qname - the id of sequence where feature was annotated
    qstart - start position of feature
    qend - end position of feature
    annot - description of the feature

    Optionally, file can store additional columns that provide additional data
    about the features.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.
    features_to_load : list
        List of features that should be loaded. Each feature will be loaded
        from the file with name pattern <feature-name>.txt. The directory to
        store files is input/annots.

    Returns:
    --------
    dict
        Dictionary where keys are names of features and values are pandas
        dataframes corresponding to these features.
    """
    features_data = {}
    for feature in features_to_load:
        ffilepath = workdir_path \
        + '/ecf-explorer/input/annots/{}.txt'.format(feature)
        fdata = pd.read_csv(ffilepath, sep=',')
        features_data[feature] = fdata

    return features_data

def select_features_to_display(features_data, reprseqs_to_display):
    """From the dataframes of annotated features, select only those rows
    that contains data for representative sequences with the numbers given in
    reprseqs_to_display list.

    Parameters:
    -----------
    features_data : dict
        Dictionary where keys are names of features and values are pandas
        dataframes corresponding to these features.
    reprseqs_to_display : list
        List of numerical identifiers of representative proteins. The input
        dataframes will be filtered so it will contain data regarding only
        those proteins.

    Returns:
    --------
    dict
        Dictionary where keys are names of features and values are pandas
        dataframes corresponding to these features. The dataframes are filtered
        and contain only data regarding only selected proteins.
    """
    selected_features_data = {}
    # if none data given: return empty dict
    if not features_data:
        return selected_features_data

    # convert list of numerical identifiers to full string identifiers
    feature_df = list(features_data.values())[0]
    rseq_id_width = int(feature_df.iloc[0].str.count(r'\d')['qname'])
    reprseq_ids = ['reprseq{:0>{idw}}'.format(rseq, idw=rseq_id_width) \
    for rseq in reprseqs_to_display]

    # select data from feature dataframes
    for fname, fdata in features_data.items():
        selected_fdata = fdata[fdata['qname'].isin(reprseq_ids)]
        selected_features_data[fname] = selected_fdata

    return selected_features_data
