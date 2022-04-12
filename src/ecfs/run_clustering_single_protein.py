import sys
import pandas as pd

if __name__ == '__main__':
    protein_df_path = sys.argv[1]
    eps             = float(sys.argv[2])
    min_samples     = float(sys.argv[3])

    sys.path.append(sys.argv[4])
    from clustering.dbscan import run_dbscan_clustering_on_single_protein_df
    from clustering.dbscan import calculate_cluster_stats

    # run dbscan on single protein df
    protein_df        = pd.read_csv(protein_df_path, sep=',')
    dbscan_params     = {'eps': eps, 'min_samples': min_samples}
    protein_df_sorted = run_dbscan_clustering_on_single_protein_df(protein_df,
                                                                  dbscan_params)

    # calculate stats for each cluster
    protein_stats = calculate_cluster_stats(protein_df_sorted)

    # save df after clustering and stats df into files
    protein_df_sorted.to_csv(protein_df_path.replace('.txt', '.clust.out'),
                             index=False)
    protein_stats.to_csv(protein_df_path.replace('.txt', '.stat.out'),
                             index=False)
