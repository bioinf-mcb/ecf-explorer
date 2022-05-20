import matplotlib.pyplot as plt

def setup_figure(protein_len, protein_id):
    """Setup figure, its axes and details (axes limits, labels, titles).

    Parameters:
    -----------
    protein_len : int
        The length of a query protein plotted.
    protein_id : string
        The identifier of a query protein plotted.

    Returns:
    --------
    fig : matplotlib.figure.Figure
        Matplotlib object representing figure (both hit alignment and features
        plots).
    ax_feat : matplotlib.axes._axes.Axes
        Matplolib object representing features subplot.
    ax_aln : matplotlib.axes._axes.Axes
        Matplolib object representing hits alignment subplot.
    """
    fig = plt.figure(figsize=(16,10))
    ax_feat = fig.add_axes([0.05, 0.85, 0.9, 0.1])
    ax_aln = fig.add_axes([0.05, 0.05, 0.9, 0.75])

    ax_feat.set_title('Features alignment to the sequence [' \
    + protein_id + ']', fontsize=16)
    ax_feat.set_xlim(0,protein_len+10)
    ax_aln.set_title('Hits alignment to the sequence [' \
    + protein_id + ']', fontsize=16)
    ax_aln.set_xlim(0,protein_len+10)
    ax_aln.set_xlabel('Query sequence [' + protein_id + ']', fontsize=14)

    return fig, ax_feat, ax_aln

def prepare_coloring(unique_labels):
    """Prepare coloring schemes for plots. Depending on the number of unique
    cluster labels select Matplotlib build-in pallete.

    Parameters:
    -----------
    unique_labels : set
        Set of unique cluster labels from protein clustering step.

    Returns:
    --------
    color2label : dict
        Dictionary mapping cluster identifier with assigned color. The color is
        represented as 3-element tuple (RGB coding).
    feature_colors : tuple
        Tuple of colors for feature visualization. Each color is represented as
        3-element tuple (RGB coding).
    """
    if len(unique_labels) <= 10:
        colors = plt.get_cmap('tab10').colors
    else:
        colors = plt.get_cmap('tab20').colors
    color2label = {}
    for k, col in zip(unique_labels, colors):
        color2label[k] = col

    feature_colors = plt.get_cmap('tab10').colors

    return color2label, feature_colors

def plot_hits(df, color2label, ax_aln, protein_id):
    """Plot hits alignment to the selected protein.

    Parameters:
    -----------
    df : pandas.DataFrame
        Pandas dataframe with results of all vs all hhblits search of
        representative proteins that was already processed with clustering.
        Filtered, so it contains only data regarding single protein to be
        plotted.
    color2label : dict
        Dictionary mapping cluster identifier with assigned color. The color is
        represented as 3-element tuple (RGB coding).
    ax_aln : matplotlib.axes._axes.Axes
        Matplolib object representing hits alignment subplot.
    protein_id : string
        The identifier of a query protein plotted.

    Returns:
    --------
    """
    for n, hit in enumerate(df.iterrows()):
        if (len(df[df['sname'] == hit[1]['sname']]) > 1):
            lstyle = "dashdot"
        else:
            lstyle = "solid"
        if protein_id == hit[1]['sname']:
            color_align = 'black'
        else:
            color_align = color2label[hit[1]['cluster']]
        ax_aln.plot([hit[1]['qstart'], hit[1]['qend']], [n,n],
        color=color_align, linestyle=lstyle)

def plot_features(selected_features_data, protein_id, ax_feat, feature_colors):
    """Plot features alignment to the selected protein.

    Parameters:
    -----------
    selected_features_data : dict
        Dictionary where keys are names of features and values are pandas
        dataframes corresponding to these features. The dataframes are filtered
        and contain only data regarding only selected proteins.
    protein_id : string
        The identifier of a query protein plotted.
    ax_feat : matplotlib.axes._axes.Axes
        Matplolib object representing features subplot.
    feature_colors : tuple
        Tuple of colors for feature visualization. Each color is represented as
        3-element tuple (RGB coding).

    Returns:
    --------
    """
    for fid, (feature_name, feature_df) in \
    enumerate(selected_features_data.items()):
        for i, hit in feature_df[feature_df['qname'] \
        == protein_id].iterrows():
            ax_feat.plot([hit['qstart'], hit['qend']], [fid, fid],
            color=feature_colors[fid], alpha=0.5, linewidth=3)
        # make space on the plot for feature not detected in a protein
        if feature_df[feature_df['qname'] == protein_id].empty:
            ax_feat.plot([0, 0], [fid, fid], color='white', alpha=1.0)
        ax_feat.text(0, fid-len(selected_features_data)*0.04,
        feature_name, fontsize=8, color=feature_colors[fid])

    return 0

def plot_alignment_selected(workdir_path, selected_reprseq_data,
    selected_features_data, save_figure=False):
    """Plot alignment of hhblits hits to a protein. Each color annotates
    different ECF. ECFs are annotated by considering each cluster as a separate
    ECF.
    Dashdot lines indicate hits to multi-hit proteins, i.e. proteins that
    are present in hit set more than once, each hit indicating different region
    of protein.
    Black lines indicate self-hits, i.e. query and subject protein is the same.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.
    selected_reprseq_data : pd.DataFrame
        Pandas dataframe with results of all vs all hhblits search of
        representative proteins that was already processed with clustering and
        filtered so only results of selected representative proteins are
        present.
    selected_features_data : dict
        Dictionary where keys are names of features and values are pandas
        dataframes corresponding to these features. The dataframes are filtered
        and contain only data regarding only selected proteins.
    save_figure : bool
        Switch on/off saving figure. The figure will be saved into plots/
        directory. The name of the figure will be reprseq<id>.png.

    Returns:
    --------
    """
    # group results by representative sequence
    reprseq_by_qname = selected_reprseq_data.groupby('qname')

    for i, df in reprseq_by_qname:
        protein_id    = df['qname'].tolist()[0]
        protein_len   = df['qlength'].tolist()[0]
        unique_labels = set(df['cluster'].tolist())

        # create figure, add axes and labels
        fig, ax_feat, ax_aln = setup_figure(protein_len, protein_id)
        # prepare coloring for alignment & features plot
        color2label, feature_colors = prepare_coloring(unique_labels)
        # plot alignments
        plot_hits(df, color2label, ax_aln, protein_id)
        # check if features are to be plotted and do it if yes
        if selected_features_data:
            plot_features(selected_features_data, protein_id, ax_feat,
                          feature_colors)
        plt.show()

        if save_figure:
            save_filename = workdir_path + '/ecf-explorer/plots/' + protein_id \
            + '.png'
            fig.savefig(save_filename)
        plt.close()
