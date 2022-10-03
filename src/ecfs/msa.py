import numpy as np
from Bio import SeqIO
from csb.bio.io import HHOutputParser

def get_qecf_hits(clustering_results, qname, cluster_id):
    """Extract hits comprising given qECF.

    Parameters:
    -----------
    clustering_results : pd.DataFrame
        Dataframe with results of all vs all hhblits search of representative
        proteins that was already processed with clustering.
    qname : str
        Identifier of the representative sequence containing qECF of interest.
    cluster_id : int
        Identifier of the cluster (i.e. qECF) within the qname representative
        sequence.

    Returns:
    --------
    list
        List of hits comprising given qECF. Each element of list is a string
        containing reprseq id of subject sequence, its start and end
        positions within query sequence. These data is separated by '_' within
        the string.
    """
    ecf_hits = clustering_results[(clustering_results['qname'] == qname) & \
               (clustering_results['cluster'] == cluster_id)]

    return ecf_hits[['sname', 'qstart', 'qend']].apply(lambda \
                   row:"_".join([str(i) for i in row]), axis=1).tolist()

def generate_msa(hhr_file, queryseq, hits, msa_params):
    """Generate mulitple sequence alignment for a representative sequence and
    a set of its hits.

    #FXME! this function should be rewritten to adhere to coding guidelines

    Parameters:
    -----------
    hhr_file : str
        Path to file hhr file containing hhblits output for query sequence.
    queryseq : str
        Sequence of the query protein.
    hits : list
        List of hits comprising given qECF. Each element of list is a string
        containing reprseq id of subject sequence, its start and end
        positions within query sequence. These data is separated by '_' within
        the string.
    msa_params : dict
        Parameters to be used to filter out hits. Available parameters:
        - ident_cut : discard hits with identity lower than the value
        - qcov_cut : discard hits with query coverage lower than the value
        - eval_cut : discard hits with e-value higher than the value

    Returns:
    --------
    list
        Multiple sequence alignment extracted from hhr file. Each element of
        the list is single sequence comprising alignment.
    list
        List of representative sequence identifiers for MSA. Each element of
        the list is identifier of one sequence. Order of identifiers
        corresponds with order of sequences in MSA list.
    """
    fasta     = [queryseq]
    hit_order = []

    for hit in HHOutputParser(alignments=True).parse_file(hhr_file):
        hit_id = '{}_{}_{}'.format(hit.id, hit.qstart, hit.qend)
        if not hit_id in hits:
            continue

        query_cov = 1.*len(hit.alignment.subject.replace('-', ''))
        if hit.identity < msa_params['ident_cut']:
            continue
        if query_cov / len(queryseq) < msa_params['qcov_cut']:
            continue
        if hit.evalue > msa_params['eval_cut']:
            continue

        hit_order.append(hit.id)
        temp  = ''
        mpos  = 0
        subject = "-"*(hit.qstart-1) + hit.alignment.subject

        # for each aa in subject
        for i in range(len(subject)):
            # no insertion at this position
            if i - hit.qstart + 1 < 0 or \
            hit.alignment.query[i - hit.qstart+ 1 ] != '-':
                while fasta[0][mpos] == '-':
                    mpos = mpos + 1
                    temp = temp + '-'
                temp = temp + str(subject[i])
                mpos = mpos + 1
            # insertion present
            else:
                if fasta[0][mpos] != '-':
                    for f in range(len(fasta)): # we need to add a gap
                        fasta[f] = fasta[f][:mpos] + "-" + fasta[f][mpos:]
                temp = temp + str(subject[i])
                mpos = mpos + 1

        fasta.append(temp)

    # fill gaps at the N terminus
    for f in range(len(fasta)):
        if len(fasta[f])<len(fasta[0]):
            fasta[f]=fasta[f]+ "-" * (len(fasta[0])-len(fasta[f]))

    return fasta, hit_order

def clean_msa(msa):
    """Remove columns in MSA where only query sequence has amino acid.

    Parameters:
    -----------
    msa : list
        Multiple sequence alignment extracted from hhr file. Each element of
        the list is single sequence comprising alignment.

    Returns:
    --------
    list
        Multiple sequence alignment extracted from hhr file. Each element of
        the list is single sequence comprising alignment. The MSA has removed
        columns where only query sequence had amino acid present.
    """

    msa_clean = []

    msa_cols = [ [ aa for aa in seq ] for seq in msa ]
    msa_arr  = np.array(msa_cols)

    # take each column (alignment position) one by one and iterate over
    # ommit first row (query row)
    for pos in msa_arr.T:
        for aa in pos[1:]:
            if aa != '-':
                msa_clean.append(pos)
                break

    # transpose back to get original array orientation
    msa_clean = np.array(msa_clean)
    msa_clean = msa_clean.T

    return msa_clean

def create_progress_bar(clustering_stats):
    """Create and print progress bar basing on number of items to process.

    Parameters:
    -----------
    clustering_stats : pd.DataFrame
        Dataframe where each row represents one cluster (qECF) and columns
        contain calculated statistics for clusters.

    Returns:
    --------
    list
        List of 100 steps that will marks progress in computation.
    """
    ecfs_n = len(clustering_stats)
    step   = int(ecfs_n / 100)
    pbar   = [ n* step for n in range(100) ] + [ecfs_n -1]
    print('[{:^100}]'.format('Building MSA for all qECFs...'))
    print('[{:^100}]'.format(str(ecfs_n) + ' qECFs to process'))

    return pbar[1:]

def show_progress_bar(pbar, n):
    """Update progress bar.

    Parameters:
    -----------
    pbar : list
        List of 100 steps that will marks progress in computation.
    n : int
        Number of already processed items.

    Returns:
    --------
    """
    if n == 0:
        print('[', end='')
    if n in pbar:
        print('=', end='')
    if n == pbar[-1]:
        print(']')

def save_qecf_as_msa(workdir_path, clustering_stats, clustering_results,
                     msa_params):
    """Save each qECF as multiple sequence alignment of its comprising hits.
    Requires repr-seqs.fa file in proper location for extraction of
    representative sequences. Requires also .hhr files from all-vs-all HHblits
    comparison of all representative sequences.

    Parameters:
    -----------
    workdir_path : string
        Path to the working directory.
    clustering_stats : pd.DataFrame
        Dataframe where each row represents one cluster (qECF) and columns
        contain calculated statistics for clusters.
    clustering_results : pd.DataFrame
        Dataframe with results of all vs all hhblits search of representative
        proteins that was already processed with clustering.
    msa_params : dict
        Parameters to be used to filter out hits. Available parameters:
        - ident_cut : discard hits with identity lower than the value
        - qcov_cut : discard hits with query coverage lower than the value
        - eval_cut : discard hits with e-value higher than the value

    Returns:
    --------
    """
    repr_seqs_filepath = workdir_path + '/ecf-explorer/input/repr-seqs.fa'
    all_by_all_hhr_dirpath = workdir_path + '/ecf-explorer/input/all-vs-all/'
    msa_output_dirpath = workdir_path + '/ecf-explorer/intermediate/qecf_msa/'

    try:
        repr_seqs = SeqIO.to_dict(SeqIO.parse(repr_seqs_filepath, 'fasta'))
    except FileNotFoundError:
        print('Please insert repr-seqs.fa file into /ecf-explorer/input/ dir')
        return False
    progress_bar = create_progress_bar(clustering_stats)

    for i, qecf in clustering_stats.iterrows():
        qname      = qecf['qname']
        cluster_id = qecf['cluster']
        queryseq   = str(repr_seqs[qname].seq)
        ecf_hits   = get_qecf_hits(clustering_results, qname, cluster_id)

        msa, hit_order = generate_msa(all_by_all_hhr_dirpath + qname + '.hhr',
                                      queryseq, ecf_hits, msa_params)
        msa_clean = clean_msa(msa)

        # add reprseq name with cluster id to list of hits
        hit_order = [qname] + hit_order
        hit_order[0] = hit_order[0] + '_' + str(cluster_id)

        if len(hit_order) >= 5:
            fmsa = open(msa_output_dirpath + 'ecf_' + qname + '_' \
                        + str(cluster_id) + '.a3m', 'w')
            for seq, hit_id in zip(msa_clean, hit_order):
                fmsa.write('>' + hit_id + '\n')
                fmsa.write(''.join(seq) + '\n')
            fmsa.close()
        show_progress_bar(progress_bar, i)
