import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN

def run_dbscan_clustering_on_single_protein_df(protein_df, dbscan_params):
    """Run DBSCAN algorithm on positions (qstart, qend) of hits to one protein
    and append to original dataframe a new column with clusters assigned to
    hits. Hits not assigned to any cluster are droped.

    Parameters:
    -----------
    protein_df : pd.DataFrame
        Pandas dataframe with parsed results from HHblits on one protein.
    dbscan_params : dict
        Dictionary of parameters used in DBSCAN (epsilon and min_samples).

    Returns:
    --------
    pd.DataFrame
        Pandas dataframe with appended column "cluster" indicating cluster
        assigned to each hit. Returned df is sorted by cluster, then qstart,
        then slength. Hits with no cluster assigned are droped. If no clusters
        were assigned returns empty dataframe.
    """
    if protein_df.empty: return None
    protein_hits = protein_df[['qstart', 'qend']]

    # data scaling
    X = protein_hits.to_numpy()
    X = StandardScaler().fit_transform(X)

    # execute DBSCAN
    eps = dbscan_params['eps']
    min_samples = dbscan_params['min_samples']
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # add to df cluster labels
    protein_df['cluster'] = labels
    # sort by cluster, start position and length (useful for plots)
    protein_df_sorted = protein_df.sort_values(by=['cluster', 'qstart',
                                                    'slength'])

    # exclude proteins without assigned clusters (considered noise)
    protein_df_sorted = protein_df_sorted[protein_df_sorted['cluster'] != -1]

    return protein_df_sorted

def calculate_cluster_stats(protein_df_sorted):
    """Calculate statistics for hit clusters created basing on positions
    (qstart, qend) of hits to the protein. Several statistics are computed:

    nb_hits - number of hits comprising cluster;
    avg_hits_overlap - average overlap of hits comprising cluster;
    ecf_span - span of cluster hits on the query protein;
    avg_hit_len_to_ecf_span - ratio of average hit length to cluster span;
    span_to_qlength - ratio of cluster span to length of query protein;

    Details regarding how those are calculated are explained in functions.

    Parameters:
    -----------
    protein_df_sorted : pd.DataFrame
        Pandas dataframe with parsed results from HHblits on one protein. The
        results are filtered and sorted by assigned cluster.

    Returns:
    --------
    pd.DataFrame / None
        Dataframe where each row represents one cluster and columns contain
        calculated statistics for a given cluster. If the input dataframe is
        empty returns empty dataframe to mark given protein as processed.
    """

    def calculate_average_hits_overlap(cluster):
        """Calculate average of overlap between all hits. Overlap is defined as
        size of intersection divided by size of union of two sets, each created
        from coords of given hit on query. If there is only one hit in a cluster
        (can happen with min_sample=1 in DBSCAN), 0.0 is returned.

        Parameters:
        -----------
        cluster : pd.DataFrame
            Pandas dataframe with parsed results from HHblits. The data
            encompass single cluster from a given protein.

        Returns:
        --------
        float
            Average overlap of hits within the cluster. 0.0 is none overlap at
            all and 1.0 is perfect overlap of all hits.
        """

        overlaps = {}

        for hid1, h1  in cluster.iterrows():
            for hid2, h2 in cluster.iterrows():
                if hid1 == hid2:
                    continue
                oid = str(hid1) + '_' + str(hid2)
                if oid not in overlaps.keys():
                    hrange1 = set(range(h1.qstart, h1.qend+1))
                    hrange2 = set(range(h2.qstart, h2.qend+1))
                    overlaps[oid] = len(hrange1.intersection(hrange2)) \
                    / len(hrange1.union(hrange2))

        overlaps_values = list(overlaps.values())
        if overlaps_values:
            return np.mean(overlaps_values)
        else:
            return 0.0

    def calculate_ECF_span(cluster):
        """Calculate the span of ECF defined as difference between largest value
        of qend position and smallest value of qstart position from all hits
        within the cluster. Returns those start/stop positions too.

        Parameters:
        -----------
        cluster : pd.DataFrame
            Pandas dataframe with parsed results from HHblits. The data
            encompass single cluster from a given protein.

        Returns:
        --------
        int
            The number of positions that comprise ECF span.
        int
            Position desginated as start of ECF.
        int
            Position desginated as end of ECF.
        """
        ecf_start = cluster.qstart.min()
        ecf_stop = cluster.qend.max()
        return ecf_stop - ecf_start + 1, ecf_start, ecf_stop

    def calculate_avg_hit_len_to_ecf_span(cluster, ecf_span):
        """Calculate ratio of average hit length to ECF span.

        Parameters:
        -----------
        cluster : pd.DataFrame
            Pandas dataframe with parsed results from HHblits. The data
            encompass single cluster from a given protein.
        ecf_span : int
            Span of ECF defined as difference between largest value of qend
            position and smallest value of qstart position from all hits within
            the cluster. Must be span corresponding to data in cluster argument.

        Returns:
        --------
        float
            Ratio of average hit length in the cluster to the span of this
            cluster.
        """

        hit_lengths = cluster['qend'] - cluster['qstart']
        avg_hit_len = hit_lengths.mean()

        return avg_hit_len / ecf_span

    def calculate_ecf_span_to_query_length(cluster, ecf_span):
        """Calculate ratio of ECF span to full query length.

        Parameters:
        -----------
        cluster : pd.DataFrame
            Pandas dataframe with parsed results from HHblits. The data
            encompass single cluster from a given protein.
        ecf_span : int
            Span of ECF defined as difference between largest value of qend
            position and smallest value of qstart position from all hits within
            the cluster. Must be span corresponding to data in cluster argument.

        Returns:
        --------
        float
            Ratio of cluster span to the full length of a query protein.
        """

        protein_len = cluster['qlength'].tolist()[0]

        return ecf_span / protein_len

    if protein_df_sorted.empty:
        return pd.DataFrame(columns=['qname', 'cluster', 'nb_hits',
                                     'avg_hits_overlap', 'ecf_span',
                                     'avg_hit_len_to_ecf_span',
                                     'span_to_qlength', 'ecf_start',
                                     'ecf_stop'])

    protein_cluster_stats = {}
    protein_id            = protein_df_sorted['qname'].tolist()[0]

    for n, (cid, cluster) in enumerate(protein_df_sorted.groupby(['cluster'])):
        avg_hits_overlap        = calculate_average_hits_overlap(cluster)
        ecf_span, ecf_start, ecf_stop = calculate_ECF_span(cluster)
        avg_hit_len_to_ecf_span = calculate_avg_hit_len_to_ecf_span(cluster,
                                                                    ecf_span)
        span_to_qlength         = calculate_ecf_span_to_query_length(cluster,
                                                                    ecf_span)

        stats = [protein_id, cid, len(cluster), avg_hits_overlap, ecf_span,
                 avg_hit_len_to_ecf_span, span_to_qlength, ecf_start, ecf_stop]
        protein_cluster_stats[n] = stats

    protein_stat_df = pd.DataFrame.from_dict(protein_cluster_stats,
                                   orient='index', columns=['qname', 'cluster',
                                   'nb_hits', 'avg_hits_overlap', 'ecf_span',
                                   'avg_hit_len_to_ecf_span',
                                   'span_to_qlength', 'ecf_start', 'ecf_stop'])

    return protein_stat_df
