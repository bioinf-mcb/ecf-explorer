import pandas as pd
import numpy as np
import sys
import pickle

if __name__ == '__main__':
	parts_dir = sys.argv[1]
	parts_out_dir = sys.argv[2]
	part_id = sys.argv[3]

	sys.path.append(sys.argv[4])
	from ecfs.ecf_dbs_compare import get_mean_number_of_hits_per_region_ecf_clust
	from ecfs.ecf_dbs_compare import get_mean_number_of_hits_per_region
	from ecfs.ecf_dbs_compare import get_continuous_regions
	from ecfs.ecf_dbs_compare import calculate_mean_region_length
	from ecfs.ecf_dbs_compare import compare_two_region_sets
	from ecfs.ecf_dbs_compare import compute_overlap_all_regions
	from ecfs.ecf_dbs_compare import update_regions_match
	from ecfs.ecf_dbs_compare import calculate_number_of_unmatched_regions
	from ecfs.ecf_dbs_compare import get_stats_reprseq

	# load dataframes & dicts
	tables = {}
	for dbid in ['ecf', 'ecod', 'pfam']:
		tables[dbid] = pd.read_csv('{}table-{}-part-{}.csv'.format(parts_dir, dbid, i))
	with open(parts_dir + 'ecf2ecf_cluster.pickle', 'rb') as handle:
		ecf2ecf_cluster = pickle.load(handle)
	with open(parts_dir + 'rseq2annot.pickle', 'rb') as handle:
		rseq2annot = pickle.load(handle)

	# get qnames set from tables
	qnames = tables['ecf']['qname'].to_list() + tables['ecod']['qname'].to_list() + tables['pfam']['qname'].to_list()
	qnames = list(set(qnames))

	stats_all = {}
	cutoffs = [0.5]
	for i, qname in enumerate(qnames):
		tables_qname = {}
		# filter table
		for tid, table in tables.items():
			tables_qname[tid] = table[table['qname'] == qname]
		# get stats
		stats = get_stats_reprseq(qname, tables_qname, cutoffs, ecf2ecf_cluster)
		# add annotation to stats
		if qname in rseq2annot.keys():
			stats.append(rseq2annot[qname])
		else:
			stats.append('Unknown')
		# add stats to dict
		stats_all[i] = stats

	columns = ['qname', 'n-ecf', 'n-ecf-unmatched',
			   'ecf-without-pfam', 'ecf-without-ecod',
			   'ecf-hit-per-reg', 'ecf-hit-per-reg-clust', 'mean-len-ecf-region',
			   'n-ecod', 'n-ecod-unmatched',
			   'ecod-without-pfam', 'ecod-without-ecf',
			   'ecod-hit-per-reg', 'mean-len-ecod-region',
			   'n-pfam', 'n-pfam-unmatched',
			   'pfam-without-ecf', 'pfam-without-ecod',
			   'pfam-hit-per-reg', 'mean-len-pfam-region',
			   'annot-len']
	for pair in ['reg-ecf-pfam-', 'reg-ecf-ecod-', 'reg-ecod-pfam-']:
		for cut in cutoffs:
			columns += [pair + str(int(cut*100))]
	columns += ['annot']

	# drop df as csv
	results = pd.DataFrame.from_dict(stats_all, orient='index', columns=columns)
	results.to_csv('{}ecf-pfam-ecod-part-{}.csv'.format(parts_out_dir, i),
				   header=False, index=False)
