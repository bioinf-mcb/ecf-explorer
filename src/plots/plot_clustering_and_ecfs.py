import matplotlib.pyplot as plt

def plot_clustering_and_ecfs_selected(selected_reprseq_data):
    """Plot alignment of hhblits hits to a protein. Each color annotates
    different ECF. ECFs are annotated by considering each cluster as a separate
    ECF.
    Dashdot lines indicate hits to multi-hit proteins, i.e. proteins that
    are present in hit set more than once, each hit indicating different region
    of protein.
    Black lines indicate self-hits, i.e. query and subject protein is the same.

    Parameters:
    -----------
    selected_reprseq_data : pd.DataFrame
        Pandas dataframe with results of all vs all hhblits search of
        representative proteins that was already processed with clustering and
        filtered so only results of selected representative proteins are
        present.

    Returns:
    --------
    """
    # group results by representative sequence
    reprseq_by_qname = selected_reprseq_data.groupby('qname')

    for i, df in reprseq_by_qname:
        protein_id    = df['qname'].tolist()[0]
        protein_len   = df['qlength'].tolist()[0]
        unique_labels = set(df['cluster'].tolist())

        # prepare coloring
        if len(unique_labels) <= 10:
            colors = plt.get_cmap('tab10').colors
        else:
            colors = plt.get_cmap('tab20').colors
        color2label = {}
        for k, col in zip(unique_labels, colors):
            color2label[k] = col

        # create figure
        fig = plt.figure(figsize=(15,30))

        # create plot and add labels
        ax2 = plt.subplot2grid((4, 4), (1, 0), colspan=4)
        ax2.set_xlim(0,protein_len+10)
        ax2.set_title('Hits alignment to the sequence', fontsize=16)
        ax2.set_xlabel('Query sequence [' + protein_id + ']', fontsize=14)

        # plot alignments
        for n, hit in enumerate(df.iterrows()):
            if (len(df[df['sname'] == hit[1]['sname']]) > 1):
                lstyle = "dashdot"
            else:
                lstyle = "solid"
            if protein_id == hit[1]['sname']:
                color_align = 'black'
            else:
                color_align = color2label[hit[1]['cluster']]
            ax2.plot([hit[1]['qstart'], hit[1]['qend']], [n,n],
            color=color_align, linestyle=lstyle)

        plt.plot()
