import sys
import re
import pandas as pd

def merge_alignment(row_data):
    """Combine parts of alignment from all rows into single string.

    Parameters:
    -----------
    row_data : list
        List containing data regarding one individual repeat. Each element of
        the list represents one parameter of the repeat. It does not contain
        data regarding repeat type that the repeat belongs to.

    Returns:
    --------
    list
        List with one element, which is full alignment of the repeat
        constructed by joining all partial alignments from different rows
        together.
    """
    rows_per_alignment = row_data.count(row_data[0]) - 1
    full_alignment = ''.join([row_data[12+i*5] for i in \
                              range(rows_per_alignment)])

    return [full_alignment]

def create_repeat_row(row, rtype_data):
    """Create complete record regarding single repeat by combining data
    regarding the repeat and the corresponding repeat type. The parameters are
    arranged in a list to correspond with columns of the pandas DataFrame that
    will be used to store all repeats data.

    Parameters:
    -----------
    row : string
        String containing data regarding one individual repeat. Parameters of
        the repeat are separated by spaces in the string. It does not contain
        data regarding repeat type that the repeat belongs to.
    rtype_data : dict
        Dictionary representing repeat type parameters. Keys are names of the
        parameters and values are values of those parameters.

    Returns:
    --------
    list
        List containing data regarding one individual repeat. Each element of
        the list represents one parameter of the repeat. It includes data
        regarding repeat type that the repeat belongs to.
    """
    row_data = [d for d in row.split(' ') if d != '']
    if len(row_data) < 12:
        print('Repeat data in file corrupted!')
        return []
    row_format = [row_data[9], rtype_data['type'], rtype_data['type-pval'],
                  rtype_data['type-len'], row_data[0]]
    row_format += [float(row_data[i]) for i in [1, 2, 3, 4]]
    qstart, qend = row_data[6].split('-')
    sstart, send = row_data[7].split('-')
    row_format += [int(i) for i in [qstart, qend, sstart, send, row_data[5]]]
    row_format += [row_data[11]]
    row_format += merge_alignment(row_data)

    return row_format

def merge_rows(repeat_rows):
    """Merge rows parsed from HHrepid that contain data regarding the same
    repeat. These do not contain data regarding repeat type that given
    repeat belongs to.

    Parameters:
    -----------
    repeat_rows : list
        List of rows parsed from HHrepid output files that contain data
        regarding individual repeats. There are multiple rows describing each
        repeat. These do not contain data regarding repeat type that given
        repeat belongs to.

    Returns:
    --------
    list
        List of rows parsed from HHrepid output files that contain data
        regarding individual repeats. There is always only one row describing
        each repeat. These do not contain data regarding repeat type that given
        repeat belongs to.
    """
    repeat_rows_merged = []
    # if input is empty (no repeats were found): return empty list
    if not repeat_rows:
        return repeat_rows_merged

    idx_all = [r.split(' ')[0] for r in repeat_rows]
    rows_per_repeat = idx_all.count(idx_all[0])
    total_repeats = int(len(repeat_rows) / rows_per_repeat)
    for i, row in enumerate(repeat_rows[:total_repeats]):
        row_full = ' '.join([repeat_rows[i+total_repeats*r] for r in \
                             range(rows_per_repeat)])
        repeat_rows_merged.append(row_full)

    return repeat_rows_merged

def parse_repeat_types(rtypes_sorted):
    """Convert raw HHrepid data grouped by repeat type into list of repeats
    that can be neatly converted into pandas DataFrame.

    Parameters:
    -----------
    rtypes_sorted : list
        List of detected repeat types in a sequence. Each type is represented
        by a list and each element of this list is a line parsed from HHrepid
        output file regarding this repeat type.

    Returns:
    --------
    list
        Data regarding individual repeats. Each repeat is represented as a list
        and each element of this list is one parameter of the repeat.
    """
    repeats_parsed = []

    for rtype in rtypes_sorted:
        rtype_data = {}
        repeat_rows = []
        for line in rtype:
            line_split = line.split(' ')
            if line_split[0] == 'Type':
                rtype_data['type'] = line_split[-1]
            if line_split[0] == 'P-value':
                rtype_data['type-pval'] = float(line_split[-1])
            elif line_split[0] == 'Length':
                rtype_data['type-len'] = int(line_split[-1])
            elif re.match(r'\D\d', line_split[0]):
                repeat_rows.append(line)
        repeat_rows = merge_rows(repeat_rows)

        repeats_parsed \
        += [create_repeat_row(row, rtype_data) for row in repeat_rows]

    return repeats_parsed

if __name__ == '__main__':
    results_out_path = sys.argv[1]

    # convert raw output from hhrepid into dataframe with repeats as rows
    rtypes_sorted = []
    with open(results_out_path, 'r') as rfile:
        # read and clean lines from file
        data = rfile.readlines()
        data = [line.strip() for line in data if line != '\n']
        current_type = []
        # group data by repeat type
        for line in data:
            if line.split(' ')[0] == 'Results':
                if current_type:
                    rtypes_sorted.append(current_type)
                current_type = ['Type ' + line.split(' ')[-1].strip(':')]
            else:
                current_type.append(line)
        rtypes_sorted.append(current_type)
    repeats_parsed = parse_repeat_types(rtypes_sorted)
    # clean corrupted entries
    repeats_parsed = [repeat for repeat in repeats_parsed if repeat]

    # save result to csv file
    columns = ['qname', 'annot', 'type-pval', 'type-len', 'rep-id',
    'prob', 'pval', 'rep-score', 'rep-score-norm', 'qstart', 'qend', 'sstart',
    'send', 'rep-len', 'rep-offset', 'rep-seq']
    repeats_df = pd.DataFrame(repeats_parsed, columns=columns)
    repeats_df.to_csv(results_out_path.replace('.out', '.csv'), index=False)
