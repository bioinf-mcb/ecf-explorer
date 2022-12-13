import os
import shutil
import pandas as pd

def setup_dir_tree(input_table_filepath):
	"""Create all directories for ECF search and analysis. They will be created
	in the same directory as the location of the input file. Created dir tree
	looks as follows:

	ecf-explorer
	|
	|-- input (copy input files here, both all-vs-all and features tables)
	|-- intermediate (df filtered by parameters and single protein dfs)
	|-- output (summary dfs with detected clusters/ECFs and their statistics)
	|-- plots (save plots here)
	|-- scripts (save bash scripts for parallel runs here)

	Parameters:
	-----------
	input_table_filepath : string
		String path with location of pandas-readable table with results of
		all vs all comparison of proteins done with HHsuite.

	Returns:
	--------
	"""
	table_dirpath = '/'.join(input_table_filepath.split('/')[:-1])

	for dir_name in ['ecf-explorer', 'ecf-explorer/input',
					 'ecf-explorer/input/annots',
					 'ecf-explorer/input/all-vs-all',
					 'ecf-explorer/input/annots/hhrepid/',
					 'ecf-explorer/intermediate',
					 'ecf-explorer/intermediate/protein_dfs',
					 'ecf-explorer/intermediate/qecf_msa',
					 'ecf-explorer/intermediate/qecf_dfs',
					 'ecf-explorer/intermediate/qecf_all_by_all',
					 'ecf-explorer/output',
					 'ecf-explorer/plots', 'ecf-explorer/scripts']:
		try:
			os.mkdir(table_dirpath + '/' + dir_name)
		except FileExistsError:
			print('/' + dir_name + '/ directory already set up')

def get_repeat_positions(data):
	"""Get all positions annotated as containing repeats in a given reprseq.

	Parameters:
	-----------
	data : pandas.DataFrame
		Dataframe containing data parsed from HHrepid output files regarding
		individual repeats. There is always only one row describing
        each repeat. This dataframe contains data reagrding only one reprseq.

	Returns:
	--------
	set
		Set of positions within reprseq that are annotated with at least one
		repeat.
	"""
	pos_with_repeat = []
	for i, repeat in data.iterrows():
		repeat_pos = list(range(repeat['qstart'], repeat['qend']+1))
		pos_with_repeat += repeat_pos

	return set(pos_with_repeat)

def check_if_repeat(hit, reprseq2repeat):
	"""Check if there is a repeat within given hit.

	Parameters:
	-----------
	hit : pandas.DataFrame
		DataFrame containing data regarding single hit to query reprseq.
	reprseq2repeat : dict
		Dictionary mapping reprseqs with the set of repeat-containing positions
		within this reprseq.

	Returns:
	--------
	bool
		Return True if there is a repeat within hit and False if there is not.
	"""
	hit_span = set(list(range(hit.qstart, hit.qend+1)))
	if hit['qname'] in reprseq2repeat.keys():
		repeats_span = reprseq2repeat[hit['qname']]
	else:
		repeats_span = set()
	common_positions = hit_span.intersection(repeats_span)

	if common_positions:
		return True
	else:
		return False

def load_and_filter_data(input_table_filepath, filters_used, filters_params):
	"""Load data from the table, filter according to parameters and return as
	a new dataframe. Additionally, save as a new dataframe in the /intermediate
	directory.

	Parameters:
	-----------
	input_table_filepath : string
		String path with location of pandas-readable table with results of
		all vs all comparison of proteins done with HHsuite.
	filters_used : list
		List of filters to be applied. To use just add a string with filter
		name, e.g. ['eval', 'self']. Available options:
		'prob' - probability value from hhblits search;
		'eval' - e-value from hhblits search;
		'self' - exclude self-hits, i.e. any hit were qname and sname are same;
		'self-full' - exclude self-hits but only if they encompass full length
		of query protein;
		'cov' - hit coverage defined as a fraction of query coveraged by hit;
		'cov-both' - hit coverage defined as both fraction of query and subject
		covered by a hit, i.e. both coverages must be met for hit to pass;
		'hit-len' - length that hit coveres on the query protein;
		'repeats' - exclude hits containing detected repeats;
	filters_params : dict
		Dictionary of values to be applied to filters used. To use add
		key-value pair to dict, e.g. {'prob_threshold':95.0}. Available options:
		'prob_threshold' <float> - minimal probability for hit to pass;
		'eval_threshold' <float> - maxminal evalue for hit to pass;
		'coverage_cutoff' <float> - minimal coverage for hit to pass, this
		parameter applies to both cov and cov-both options;
		'min_hit_len' <int> - minimal length that hit coveres on the query
		protein;
		'repeat_pval' <float> - maximum p-value for repeat to be considered;

	Returns:
	--------
	pandas.DataFrame
		Pandas dataframe filtered according to the selected parameters. The
		dataframe is also saved in the /intermediate directory.
	"""
	###  setup paths, copy input and load data
	input_table_dirpath  = '/'.join(input_table_filepath.split('/')[:-1]) \
	+ '/ecf-explorer/input/'
	output_table_dirpath = '/'.join(input_table_filepath.split('/')[:-1]) \
	+ '/ecf-explorer/intermediate/'
	input_filename       = input_table_filepath.split('/')[-1]

	shutil.copyfile(input_table_filepath, input_table_dirpath + input_filename)
	hhr_table = pd.read_csv(input_table_dirpath + input_filename, sep=',')

	### filter data according to previously set parameters
	# by probability
	if 'prob' in filters_used:
		hhr_table = hhr_table[hhr_table['prob'] >= \
		filters_params['prob_threshold']]
	# by evalue
	if 'eval' in filters_used:
		hhr_table = hhr_table[hhr_table['eval'] <= \
		filters_params['eval_threshold']]

	# by coverage
	if 'cov' in filters_used:
		hhr_table = hhr_table.assign(cov = lambda x: \
		(x.qend - x.qstart + 1) / x.qlength)
		hhr_table = hhr_table[hhr_table['cov'] \
		>= filters_params['coverage_cutoff']]
	# by-directional coverage, qcov and scov
	if 'cov-both' in filters_used:
		hhr_table = hhr_table.assign(qcov = lambda x: \
		(x.qend - x.qstart + 1) / x.qlength)
		hhr_table = hhr_table.assign(scov = lambda x: \
		(x.send - x.sstart + 1) / x.slength)
		hhr_table = hhr_table[hhr_table['qcov'] \
		>= filters_params['coverage_cutoff']]
		hhr_table = hhr_table[hhr_table['scov'] \
		>= filters_params['coverage_cutoff']]
	# by len of query covered by hit
	if 'hit-len' in filters_used:
		hhr_table = hhr_table.assign(hlen = lambda x: (x.qend - x.qstart + 1))
		hhr_table = hhr_table[hhr_table['hlen'] \
		>= filters_params['min_hit_len']]

	# exclude self-hits
	if 'self' in filters_used:
		hhr_table = hhr_table[hhr_table['sname'] != hhr_table['qname']]
	# exclude self-hits but only if they encompass full protein
	if 'self-full' in filters_used:
		hhr_self          = hhr_table[hhr_table['sname'] == hhr_table['qname']]
		hhr_self_not_full = hhr_self[(hhr_self['qend'] - hhr_self['qstart'] + 1)
		!= hhr_self['qlength']]
		hhr_table         = hhr_table[hhr_table['sname'] != hhr_table['qname']]
		hhr_table         = pd.concat([hhr_self_not_full, hhr_table])
	# exclude hits containing repeats
	if 'repeats' in filters_used:
		repeats_filepath = '/'.join(input_table_filepath.split('/')[:-1]) \
		+ '/ecf-explorer/input/annots/hhrepid.txt'
		repeats_table = pd.read_csv(repeats_filepath, sep=',')

		repeats_table = repeats_table[repeats_table['pval'] \
		<= filters_params['repeat_pval']]

		reprseq2repeat = {qname:get_repeat_positions(data) for qname, data in \
						  repeats_table.groupby('qname')}
		hhr_table['repeat'] = \
		hhr_table.apply(lambda hit: check_if_repeat(hit, reprseq2repeat),
						axis=1)
		hhr_table = hhr_table[hhr_table['repeat'] == False]

	### save to intermediate dir
	hhr_table.to_csv(output_table_dirpath + input_filename.split('.')[0] \
	+ '_filtered.txt', sep=',', index=False)

	return hhr_table
