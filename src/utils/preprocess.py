import os
import shutil
import pandas as pd

def setup_dir_tree(input_table_filepath):
	"""Create all directories for ECF search and analysis. They will be created
	in the same directory as the location of the input file. Created dir tree
	looks as follows:

	ecf-explorer
	|
	|-- input (copy input files here)
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
					 'ecf-explorer/intermediate',
					 'ecf-explorer/intermediate/protein_dfs',
					 'ecf-explorer/output',
					 'ecf-explorer/plots', 'ecf-explorer/scripts']:
		try:
			os.mkdir(table_dirpath + '/' + dir_name)
		except FileExistsError:
			print('/' + dir_name + '/ directory already set up')

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
	filters_params : dict
		Dictionary of values to be applied to filters used. To use add
		key-value pair to dict, e.g. {'prob_threshold':95.0}. Available options:
		'prob_threshold' <float> - minimal probability for hit to pass;
		'eval_threshold' <float> - minimal evalue for hit to pass;
		'coverage_cutoff' <float> - minimal coverage for hit to pass, this
		parameter applies to both cov and cov-both options;
		'min_hit_len' <int> - minimal length that hit coveres on the query
		protein;

	Returns:
	--------
	pd.DataFrame
		Pandas dataframe filtered according to the selected parameters. The
		dataframe is also saved in the /intermediate direcotry.
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

	### save to intermediate dir
	hhr_table.to_csv(output_table_dirpath + input_filename.split('.')[0] \
	+ '_filtered.txt', sep=',', index=False)

	return hhr_table
