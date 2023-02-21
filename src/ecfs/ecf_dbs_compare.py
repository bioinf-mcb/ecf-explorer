import pandas as pd
import numpy as np

def get_mean_number_of_hits_per_region_ecf_clust(regions, table, ecf_clusters):
	"""Calculate mean number of hits comprising detected ECF region.
	Assume ECFs from the same cluster count as one hit."""
	if table.empty:
		return np.NaN

	hit_per_region = []
	processed_hits = []
	for region in regions:
		reg_hits = []
		reg_pos = set(region)
		for hid, hit in table.iterrows():
			if hid in processed_hits:
				continue
			else:
				hit_pos = set(range(hit['qstart'], hit['qend']+1))
			if reg_pos.intersection(hit_pos):
				processed_hits.append(hid)

				reg_hits.append(ecf_clusters[hit['sname']])
		hit_per_region.append(len(set(reg_hits)))

	if hit_per_region:
		mean_hit_per_region = np.mean(hit_per_region)
	else:
		mean_hit_per_region = np.NaN

	return mean_hit_per_region

def get_mean_number_of_hits_per_region(regions, table):
	"""Calculate mean number of hits comprising detected region."""
	if table.empty:
		return np.NaN
	else:
		return len(table) / len(regions)

def get_continuous_regions(table):
	"""Define continuous regions of hits(ECF, Pfam, ECOD) on reprseq."""
	if table.empty:
		return []

	positions = []
	regions = []
	for hid, hit in table.iterrows():
		positions += list(range(hit['qstart'], hit['qend']+1))
	positions = sorted(list(set(positions)))

	region = []
	for pos in range(len(positions)-1):
		if positions[pos]+1 == positions[pos+1]:
			region.append(positions[pos])
		else:
			region.append(positions[pos])
			regions.append(region)
			region = []
	regions.append(region)
	regions = [r for r in regions if r]

	return regions

def calculate_mean_region_length(regions):
	if regions:
		return np.mean([len(r) for r in regions])
	else:
		return np.NaN

def compare_two_region_sets(regions1, regions2, cutoff):
	"""Compare two sets of regions."""

	reg1_to_match = [0 for i in range(len(regions1))]
	reg2_to_match = [0 for i in range(len(regions2))]

	for rid1, reg1 in enumerate(regions1):
		for rid2, reg2 in enumerate(regions2):
			reg1_set = set(reg1)
			reg2_set = set(reg2)
			regs_union = reg1_set.union(reg2_set)
			regs_intersect = reg1_set.intersection(reg2_set)
			coverage = len(regs_intersect) / len(regs_union)
			if coverage >= cutoff:
				reg1_to_match[rid1] = 1
				reg2_to_match[rid2] = 1

	max_possible_match = len(regions1) + len(regions2)
	total_matched = sum(reg1_to_match) + sum(reg2_to_match)
	match_ratio = total_matched / max_possible_match
	match_status = [reg1_to_match, reg2_to_match, reg1_to_match.count(0), reg2_to_match.count(0)]

	return match_ratio, match_status

def compute_overlap_all_regions(regions_all, qlen):
	if qlen == 0.0:
		return 0

	regions_concat = []
	for reg_table in regions_all:
		for reg in reg_table:
			regions_concat += reg
	regions_concat = set(regions_concat)

	return round(len(regions_concat) / qlen, 3)

def update_regions_match(regions_matched, match_status, update_with_db):
	"""Update status of matched regions for a given db."""
	regions_matched_updated = []
	for reg_stat, reg_upd in zip(regions_matched, match_status):
		upd_region = reg_stat
		if reg_upd == 1:
			upd_region[update_with_db] = reg_upd
		regions_matched_updated.append(upd_region)

	return regions_matched_updated

def calculate_number_of_unmatched_regions(regions_matched):
	regions_unmatched = {}
	for db in ['ecf', 'ecod', 'pfam']:
		if regions_matched[db]:
			regions_unmatched[db] = regions_matched[db].count(0)
		else:
			regions_unmatched[db] = np.NaN

	return regions_unmatched

def calculate_venn_stats(regions_matched):
	venn_stats = {'ecf':0, 'ecf+pfam':0, 'ecf+ecod':0, 'ecf+pfam+ecod':0,
				  'pfam':0, 'pfam+ecod':0,
				  'ecod':0}

	for reg in regions_matched['ecf']:
		if (reg['pfam'] == 0 and reg['ecod'] == 0):
			venn_stats['ecf'] += 1
		elif (reg['pfam'] > 0 and reg['ecod'] == 0):
			venn_stats['ecf+pfam'] += 1
		elif (reg['pfam'] == 0 and reg['ecod'] > 0):
			venn_stats['ecf+ecod'] += 1
		elif (reg['pfam'] > 0 and reg['ecod'] > 0):
			venn_stats['ecf+pfam+ecod'] += 1

	for reg in regions_matched['pfam']:
		if (reg['ecf'] == 0 and reg['ecod'] == 0):
			venn_stats['pfam'] += 1
		elif (reg['ecf'] == 0 and reg['ecod'] > 0):
			venn_stats['pfam+ecod'] += 1

	for reg in regions_matched['ecod']:
		if (reg['ecf'] == 0 and reg['pfam'] == 0):
			venn_stats['ecod'] += 1

	return venn_stats

def get_stats_reprseq(qname, tables, cutoffs, ecf2ecf_cluster):

	# get regions
	ecf_regions = get_continuous_regions(tables['ecf'])
	ecf_hit_per_reg = get_mean_number_of_hits_per_region(ecf_regions, tables['ecf'])
	ecf_hit_per_reg_clust = get_mean_number_of_hits_per_region_ecf_clust(ecf_regions, tables['ecf'],
																		 ecf2ecf_cluster)
	ecod_regions = get_continuous_regions(tables['ecod'])
	ecod_hit_per_reg = get_mean_number_of_hits_per_region(ecod_regions, tables['ecod'])
	pfam_regions = get_continuous_regions(tables['pfam'])
	pfam_hit_per_reg = get_mean_number_of_hits_per_region(pfam_regions, tables['pfam'])
	regions_all = [ecf_regions, ecod_regions, pfam_regions]

	# get number of and mean len of regions
	n_ecf = len(ecf_regions)
	n_ecod = len(ecod_regions)
	n_pfam = len(pfam_regions)
	mean_len_ecf_region = calculate_mean_region_length(ecf_regions)
	mean_len_ecod_region = calculate_mean_region_length(ecod_regions)
	mean_len_pfam_region = calculate_mean_region_length(pfam_regions)

	# compare all regions to each other
	if ecf_regions:
		qlen = tables['ecf'].iloc[0]['qlength']
	elif ecod_regions:
		qlen = tables['ecod'].iloc[0]['qlength']
	elif pfam_regions:
		qlen = tables['pfam'].iloc[0]['qlength']
	else:
		qlen = 0.0

	regions_comparisons = {}
	regions_matched = {'ecf':[{'pfam':0, 'ecod':0} for i in range(len(ecf_regions))],
					   'ecod':[{'pfam':0, 'ecf':0} for i in range(len(ecod_regions))],
					   'pfam':[{'ecf':0, 'ecod':0} for i in range(len(pfam_regions))]}
	region_pairs_matching = {}
	for cutoff in cutoffs:
		if (ecf_regions and pfam_regions):
			regions_comparisons['reg_ecf_pfam-{}'.format(cutoff)], match_status = compare_two_region_sets(ecf_regions,
																							pfam_regions,
																							cutoff=cutoff)
			### !!! WARNING this works only for single cutoff threshold
			### !!! WARNING if more selected the data will be overwritten and reflect last cutoff
			regions_matched['ecf'] = update_regions_match(regions_matched['ecf'], match_status[0], 'pfam')
			regions_matched['pfam'] = update_regions_match(regions_matched['pfam'], match_status[1], 'ecf')
			region_pairs_matching['ecf-without-pfam'] = match_status[2]
			region_pairs_matching['pfam-without-ecf'] = match_status[3]
		else:
			regions_comparisons['reg_ecf_pfam-{}'.format(cutoff)] = np.NaN
			if ecf_regions:
				region_pairs_matching['ecf-without-pfam'] = n_ecf
				region_pairs_matching['pfam-without-ecf'] = np.NaN
			elif pfam_regions:
				region_pairs_matching['ecf-without-pfam'] = np.NaN
				region_pairs_matching['pfam-without-ecf'] = n_pfam
			else:
				region_pairs_matching['ecf-without-pfam'] = np.NaN
				region_pairs_matching['pfam-without-ecf'] = np.NaN

		if (ecf_regions and ecod_regions):
			regions_comparisons['reg_ecf_ecod-{}'.format(cutoff)], match_status = compare_two_region_sets(ecf_regions,
																							ecod_regions,
																							cutoff=cutoff)
			regions_matched['ecf'] = update_regions_match(regions_matched['ecf'], match_status[0], 'ecod')
			regions_matched['ecod'] = update_regions_match(regions_matched['ecod'], match_status[1], 'ecf')
			region_pairs_matching['ecf-without-ecod'] = match_status[2]
			region_pairs_matching['ecod-without-ecf'] = match_status[3]
		else:
			regions_comparisons['reg_ecf_ecod-{}'.format(cutoff)] = np.NaN
			if ecf_regions:
				region_pairs_matching['ecf-without-ecod'] = n_ecf
				region_pairs_matching['ecod-without-ecf'] = np.NaN
			elif ecod_regions:
				region_pairs_matching['ecf-without-ecod'] = np.NaN
				region_pairs_matching['ecod-without-ecf'] = n_ecod
			else:
				region_pairs_matching['ecf-without-ecod'] = np.NaN
				region_pairs_matching['ecod-without-ecf'] = np.NaN

		if (ecod_regions and pfam_regions):
			regions_comparisons['reg_ecod_pfam-{}'.format(cutoff)], match_status = compare_two_region_sets(ecod_regions,
																							 pfam_regions,
																							 cutoff=cutoff)
			regions_matched['ecod'] = update_regions_match(regions_matched['ecod'], match_status[0], 'pfam')
			regions_matched['pfam'] = update_regions_match(regions_matched['pfam'], match_status[1], 'ecod')
			region_pairs_matching['ecod-without-pfam'] = match_status[2]
			region_pairs_matching['pfam-without-ecod'] = match_status[3]
		else:
			regions_comparisons['reg_ecod_pfam-{}'.format(cutoff)] = np.NaN
			if ecod_regions:
				region_pairs_matching['ecod-without-pfam'] = n_ecod
				region_pairs_matching['pfam-without-ecod'] = np.NaN
			elif pfam_regions:
				region_pairs_matching['ecod-without-pfam'] = np.NaN
				region_pairs_matching['pfam-without-ecod'] = n_pfam
			else:
				region_pairs_matching['ecod-without-pfam'] = np.NaN
				region_pairs_matching['pfam-without-ecod'] = np.NaN

	# add % of qlen occupied by any db hit
	annot_region_len = compute_overlap_all_regions(regions_all, qlen)
	venn_stats = calculate_venn_stats(regions_matched)

	# calculate number of unmatched regions for each db - no longer in use
	# regions_unmatched = calculate_number_of_unmatched_regions(regions_matched)

	stats = [qname, n_ecf,
			 region_pairs_matching['ecf-without-pfam'], region_pairs_matching['ecf-without-ecod'],
			 ecf_hit_per_reg, ecf_hit_per_reg_clust, mean_len_ecf_region,
			 n_ecod,
			 region_pairs_matching['ecod-without-pfam'], region_pairs_matching['ecod-without-ecf'],
			 ecod_hit_per_reg, mean_len_ecod_region,
			 n_pfam,
			 region_pairs_matching['pfam-without-ecf'], region_pairs_matching['pfam-without-ecod'],
			 pfam_hit_per_reg, mean_len_pfam_region,
			 venn_stats['ecf'], venn_stats['ecf+pfam'], venn_stats['ecf+ecod'], venn_stats['ecf+pfam+ecod'],
			 venn_stats['pfam'], venn_stats['pfam+ecod'],
			 venn_stats['ecod'],
			 annot_region_len]
	for pair in ['reg_ecf_pfam-', 'reg_ecf_ecod-', 'reg_ecod_pfam-']:
		for cut in cutoffs:
			stats += [round(regions_comparisons[pair + str(cut)],3)]
	return stats
