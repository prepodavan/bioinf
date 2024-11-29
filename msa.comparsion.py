ascii_zero = 48
ascii_last = 126
utf_shift = 128


def score_chr(score) -> str:
    """
    Maps columns score to single character (stockholm annotations allows only single char)
    :return: printable ascii or utf-8
    :raises ValueError: if unexpected type, got non-single-char-str, got negative score or bigger than printable utf-8
    """
    if score is True:
        return '1'
    if score is False or score is None:
        return '-'
    if isinstance(score, str):
        assert len(score) == 1
        return score
    if isinstance(score, int) or isinstance(score, float):
        assert 0 <= score <= 0x10ffff
        score = int(score)
        # add code of ascii zero
        # to print ascii digits and letters
        if score + ascii_zero < ascii_last:
            score += ascii_zero
        else:
            score += utf_shift
        return chr(score)

    raise ValueError(f'unsupported type for stockholm score annotation: {score}')


aa_classes = {
    'A': 'F',
    'V': 'F',
    'L': 'F',
    'I': 'F',
    'J': 'F',
    'M': 'F',
    'P': 'F',
    'F': 'R',
    'Y': 'R',
    'W': 'R',
    'G': 'Z',
    'S': 'Z',
    'T': 'Z',
    'C': 'Z',
    'Q': 'Z',
    'N': 'Z',
    'E': 'N',
    'D': 'N',
    'K': 'P',
    'R': 'P',
    'H': 'P',
}


for residue, cls in list(aa_classes.items()):
    aa_classes[residue.lower()] = cls


def are_residues_eq(x, y) -> bool:
    """
    :return: True if X and Y equals
    """
    return x in aa_classes and x.upper() == y.upper()


def are_residues_similar(x, y) -> bool:
    """
    :return: True if X and Y have similar properties
    """
    x = aa_classes.get(x)
    y = aa_classes.get(y)
    return x is not None and x == y


class SubsMatCmp:
    """
    Scores residues during pairwise based on Bio.Align.substitution_matrices
    Example: SubsMatCmp(substitution_matrices.load("BLOSUM62")).cmp('A', 'W') == -3.0
    """
    def __init__(self, matrix):
        self.mx = matrix

    def cmp(self, x, y) -> int:
        return int(self.mx[x.upper()][y.upper()])


def fill_indices(record) -> int:
    """
    For every residue fills letter annotation with original sequence index
    :param record: SeqRecord, with or without gaps
    :return: number of last residue
    """
    # should not be printed into file and just store indices
    record.letter_annotations['fbb_indices'] = [0] * len(record.seq)
    original_index = 0
    for seq_index in range(len(record.seq)):
        if record.seq[seq_index] == '-':
            continue
        original_index += 1
        record.letter_annotations['fbb_indices'][seq_index] = original_index
    return original_index


def fill_maparray(msa1, msa2, row) -> int:
    """
    Fills veralign maparray into letter_notations.
    Maparray contains columns indexes of same alpha-C number
    :return: number of found matched residues
    """
    matches = 0
    f, s = 0, 0
    msa1[row].letter_annotations['fbb_maparray'] = [0] * len(msa1[row])
    msa2[row].letter_annotations['fbb_maparray'] = [0] * len(msa2[row])
    while f <= len(msa1[row])-1 and s <= len(msa2[row])-1:
        res1 = msa1[row].letter_annotations['fbb_indices'][f]
        res2 = msa2[row].letter_annotations['fbb_indices'][s]
        if res1 != 0 and res1 == res2:
            matches += 1
            msa1[row].letter_annotations['fbb_maparray'][f] = s+1
            msa2[row].letter_annotations['fbb_maparray'][s] = f+1
        if f == len(msa1[row])-1 or s == len(msa2[row])-1:
            break
        if res1 == 0:
            f += 1
        if res2 == 0:
            s += 1
        if res1 != 0 and res2 != 0:
            f += 1
            s += 1
    return matches


def fill_test_with_truth(test_record, truth_record):
    """
    Fills test_rec with truth_rec into letter_annotations
    :return: mutated test_rec
    """
    if len(truth_record.seq) >= len(test_record.seq):
        test_record.letter_annotations['FBB_TR'] = str(truth_record.seq)[:len(test_record.seq)]
    elif len(truth_record.seq) < len(test_record.seq):
        residue = len(test_record.seq) - len(truth_record.seq)
        seq = str(truth_record.seq)
        seq += '_' * residue
        test_record.letter_annotations['FBB_TR'] = seq
    return test_record


def column_wise(
        test_msa,
        truth_msa,
        col_annotation: str,
        aa_cmp,
        msa_annotation_prefix: str = '',
        scorech=score_chr,
    ):
    """
    Compares corresponding columns with each other and writes bool score as column_annotations to test_msa.
    :param test_msa: msa where to save scores
    :param truth_msa: msa against which score will be calculated
    :param col_annotation: name of annotation to save scores
    :param msa_annotation_prefix: if given then writes count and percentage of scores (which tests as True) to annotations of test_sma
    :param aa_cmp: comparator of residues, takes two aa and returns bool
    :return: tuple of count and percentage of scores which tests as True
    """

    assert len(test_msa) == len(truth_msa)

    matches = 0
    alnlen = test_msa.get_alignment_length()
    test_msa.column_annotations[col_annotation] = [scorech(False)] * alnlen
    for col in range(alnlen):
        got = None
        equals = True
        for row in range(len(test_msa)):
            tst = test_msa[row].seq[col].upper()
            trh_col = test_msa[row].letter_annotations['fbb_maparray'][col]
            if got is None:
                got = trh_col
            if got != trh_col:
                equals = False
                break
            if trh_col == 0:
                equals = False
                break
            
            # substract 1 since maparray starts from 1
            trh = truth_msa[row].seq[trh_col-1].upper()
            if not aa_cmp(tst, trh):
                equals = False
                break
        
        if equals:
            matches += 1
        test_msa.column_annotations[col_annotation][col] = scorech(equals)

    pcnt = matches * 100 / alnlen
    test_msa.column_annotations[col_annotation] = ''.join(test_msa.column_annotations[col_annotation])
    if len(msa_annotation_prefix) != 0:
        test_msa.annotations[f'{msa_annotation_prefix}_cw_cnt'] = matches
        test_msa.annotations[f'{msa_annotation_prefix}_cw_pcnt'] = pcnt
    return matches, pcnt


def pairwise(
        test_msa,
        truth_msa,
        col_annotation: str,
        aa_cmp,
        threshold_index: str,
        threshold = lambda t: t > 0,
        msa_annotation_prefix: str= '',
        scorech=score_chr,
    ):
    """
    Compares corresponding cells of corresponding columns with each other and writes scores as column_annotations to test_msa.
    :param test_msa: msa where to save scores
    :param truth_msa: msa against which score will be calculated
    :param col_annotation: name of annotation to save scores
    :param msa_annotation_prefix: if given then writes sum of scores (which tests as True) to annotations of test_sma
    :param aa_cmp: comparator of residues, takes two aa and returns integer score for pair
    :param threshold: function that takes column sum of pairs and returns True if threshold matched
    :param threshold_index: where to write conservativeness
    :return: tuple (sum of column scores, amount of pairs, amount of pairs with positive score)
    """

    assert len(test_msa) == len(truth_msa)

    sp = 0
    pairs = 0
    positive_pairs = 0

    alnlen = test_msa.get_alignment_length()
    test_msa.column_annotations[col_annotation] = [scorech(False)] * alnlen
    test_msa.column_annotations[threshold_index] = [scorech(False)] * alnlen
    for col in range(alnlen):
        colsum = 0
        got = None
        equals = True
        for row in range(len(test_msa)):
            pairs += 1
            tst = test_msa[row].seq[col].upper()
            trh_col = test_msa[row].letter_annotations['fbb_maparray'][col]
            if trh_col == 0:
                continue
            if got is None:
                got = trh_col
            if got != trh_col:
                equals = False
                break

            
            # substract 1 since maparray starts from 1
            trh = truth_msa[row].seq[trh_col-1].upper()
            score = int(aa_cmp(tst, trh))
            if score > 0:
                positive_pairs += 1
            colsum += score

        if not equals:
            continue

        sp += colsum
        if threshold(colsum):
            test_msa.column_annotations[threshold_index][col] = scorech(True)
        if colsum >= 0:
            test_msa.column_annotations[col_annotation][col] = scorech(colsum)

    test_msa.column_annotations[col_annotation] = ''.join(test_msa.column_annotations[col_annotation])
    test_msa.column_annotations[threshold_index] = ''.join(test_msa.column_annotations[threshold_index])
    if len(msa_annotation_prefix) != 0:
        test_msa.annotations[f'{msa_annotation_prefix}_pw_sum'] = sp
        test_msa.annotations[f'{msa_annotation_prefix}_pw_amount'] = pairs
        test_msa.annotations[f'{msa_annotation_prefix}_pw_matched'] = positive_pairs
        test_msa.annotations[f'{msa_annotation_prefix}_sp_score'] = positive_pairs * 100 / pairs
    return sp, pairs, positive_pairs


def mark_case_by_conservateveness(msa, index_name: str):
    """
    Changes cases of all sequences in msa corresponding to boolean index
    """
    from Bio.Seq import Seq


    for aln in msa:
        seq = list(aln.seq)
        for i in range(len(seq)):
            if msa.column_annotations[index_name][i] != score_chr(True):
                seq[i] = seq[i].lower()
            else:
                seq[i] = seq[i].upper()
        aln.seq = Seq(''.join(seq))
    return msa


def find_conservative_ranges(msa, index_name: str):
    """
    Finds all ranges of conservative columns.
    :param index_name: name of boolean index which determines wheteher column is conservative
    :return: list of tuples (start_test, end_test, start_truth, end_truth) of each conservative.
    """
    # could decide column is conservative
    # even if there is gapes.
    # need to find first row without gap
    # to find out position in truth alignment
    def get_maparray(column):
        for row in range(len(msa)):
            maparray = msa[row].letter_annotations['fbb_maparray']
            if len(index) > len(maparray) or maparray[column] == 0:
                continue
            return maparray[column]

    ranges = []
    index = msa.column_annotations[index_name]
    flag = False
    start = 0
    for col in range(len(index)):
        if not flag and index[col] == score_chr(True):
            flag = True
            start = col
        elif flag and index[col] == score_chr(False):
            flag = False
            trh = get_maparray(col-1)
            start_trh = get_maparray(start)
            # should not be marked as conservative in any case
            assert trh is not None
            assert start_trh is not None
            # maparray positions starts with 1 so no need to add
            ranges.append((start+1, col, start_trh, trh))

    return ranges


def load_substitution_matrices() -> list:
    """
    Returns list of tuples (matix_name, matrix) loaded from from Bio.Align.substitution_matrices.
    Matrices having not all of residues will be filtered out
    """
    from Bio.Align import substitution_matrices as sm


    out = [(name, sm.load(name)) for name in sm.load()]
    residues = set([r.upper() for r in aa_classes.keys()])
    def check(item):
        for aa in residues:
            if aa not in item[1].alphabet and aa != 'J':
                return False
        return True

    return sorted(filter(check, out), key=lambda x: x[0])


def compare_alignments(query, subject, threshold: int = 0, mark: str = None, sub_mxs=[]):
    """
    :return: tuple (mutated_query, conservative_ranges or None if mark index not given)
    """
    query.sort(key=lambda x: x.id)
    subject.sort(key=lambda x: x.id)

    for record in query:
        fill_indices(record)
    for record in subject:
        fill_indices(record)
    for i in range(min(len(query), len(subject))):
        fill_maparray(query, subject, i)
        fill_test_with_truth(query[i], subject[i])

    column_wise(query, subject, 'FBB_EQ', are_residues_eq, msa_annotation_prefix='FBB_EQ')
    column_wise(query, subject, 'FBB_SM', are_residues_similar, msa_annotation_prefix='FBB_SM')

    pairwise(
        query,
        subject,
        'FBB_SP_EQ',
        lambda x, y: int(are_residues_eq(x, y)),
        'FBB_SP_EQ_TH',
        msa_annotation_prefix='FBB_SP_EQ',
    )
    pairwise(
        query,
        subject,
        'FBB_SP_SM',
        lambda x, y: int(are_residues_similar(x, y)),
        'FBB_SP_SM_TH',
        lambda x: x >= int(threshold),
        msa_annotation_prefix='FBB_SP_SM',
    )

    for sm in sub_mxs:
        comparator = SubsMatCmp(sm[1])
        pairwise(
            query,
            subject,
            f'FBB_SP_{sm[0]}',
            comparator.cmp,
            f'FBB_SP_{sm[0]}_TH',
            msa_annotation_prefix=f'FBB_SP_{sm[0]}',
        )

    if not isinstance(mark, str) or not len(mark) != 0:
        return query, None

    mark_case_by_conservateveness(query, mark)
    ranges = find_conservative_ranges(query, mark)

    return query, ranges


def compare_alignments_files(
        query=None,
        subject=None,
        output=None,
        query_format=None,
        subject_format=None,
        output_stats=None,
        output_matrices=None,
        **kwargs,
    ):
    import io
    from Bio import AlignIO
    from Bio.AlignIO.StockholmIO import StockholmWriter


    qry = AlignIO.read(query, query_format)
    sbj = AlignIO.read(subject, subject_format)

    if output_matrices and 'sub_mxs' not in kwargs:
        kwargs['sub_mxs'] = load_substitution_matrices()
    qry, ranges = compare_alignments(qry, sbj, **kwargs)
    if ranges:
        print("Position in query,Position in subject")
        for rng in ranges:
            print(f'{rng[0]}-{rng[1]},{rng[2]}-{rng[3]}')

    # biopython 1.83 stockholm writer lacks
    # support of writing (not reading though)
    # file-level GF annotatioins.
    # biopython required that column_annotations
    # aligned by length with alignment length,
    # and letter_annotations aligned by length
    # with length of annotated sequence.
    # also, by format spec, file starts with version info,
    # and ends with '//\n'.
    # so we end up with only option - catch all writes
    # into buffer, insert stats into GF annotations
    # between header (with version and etc) and sequences
    # and only then write it to the file.
    buf = io.StringIO()
    writer = StockholmWriter(buf)
    for record in qry:
        for annotation_name in record.letter_annotations.keys():
            if annotation_name.startswith('FBB_'):
                writer.pfam_gr_mapping[annotation_name] = annotation_name
    for annotation_name in qry.column_annotations.keys():
        if annotation_name.startswith('FBB_'):
            writer.pfam_gc_mapping[annotation_name] = annotation_name
    writer.write_alignment(qry)

    out = buf.getvalue()
    if output_stats:
        stats = []
        for key, value in sorted(qry.annotations.items()):
            stats.append(f'#=GF {key} {value}\n')
        stats = ''.join(stats)

        view = buf.getvalue()
        fst_nl = view.find('\n')
        snd_nl = view.find('\n', fst_nl+1)
        out = view[:snd_nl] + "\n" + stats + view[snd_nl+1:]

    with open(output, 'w') as f:
        f.write(out)

    return qry, sbj


if __name__ == '__main__' and '__file__' in globals():
    import argparse


    parser = argparse.ArgumentParser(
        description='This program compares two MSA and writes statistics to stockholm output.',
    )
    parser.add_argument(
        "--query",
        "-r",
        required=True,
        help="alignment to test against \"the truth\"",
    )
    parser.add_argument(
        "--subject",
        "-t",
        required=True,
        help="alignment to compare with",
    )
    parser.add_argument(
        "--query-format",
        "-q",
        dest='query_format',
        required=True,
        help="format of query. required by Bio.AlignIO",
    )
    parser.add_argument(
        "--subject-format",
        "-s",
        dest='subject_format',
        required=True,
        help="format of subject. required by Bio.AlignIO",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="path to file to write output",
    )
    parser.add_argument(
        "--threshold",
        "-d",
        required=True,
        help="threshold for pairwise conservative detection"
    )
    parser.add_argument(
        "--mark-conservative-by-case",
        "-m",
        dest='mark',
        help="index name which will store sequences formated with lower case," +
            "if not conservative column, and upper case otherwise." +
            "stores nothing if not given",
    )
    parser.add_argument(
        "--stats",
        "-a",
        default=False,
        action='store_true',
        dest='output_stats',
        help="whether should output statistics"
    )
    parser.add_argument(
        "--matrices",
        "-i",
        default=False,
        action='store_true',
        dest='output_matrices',
        help="whether should build indices based on substitution matrices"
    )

    args = parser.parse_args()
    compare_alignments_files(**vars(args))

