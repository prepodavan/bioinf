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


def chr_to_score(c: str) -> int:
    """
    Maps written score character to original number
    """
    assert len(c) == 1
    score = ord(c)
    if score - ascii_zero < ascii_last:
        return score - ascii_zero
    return score - utf_shift


aa_classes = {
    'A': 'F',
    'V': 'F',
    'L': 'F',
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


def are_residues_similar(x, y):
    """
    :return: True if X and Y have same properties
    """
    x = aa_classes.get(x)
    y = aa_classes.get(y)
    return x is not None and x == y


class ConservativenessThresholdMatcher(object):
    def __init__(self, threshold, index_name):
        self.threshold = threshold
        self.index_name = index_name

    def cmp(self, test_msa, _, column: int) -> bool:
        score = test_msa.column_annotations[self.index_name][column]
        return chr_to_score(score) >= self.threshold


def column_wise(test_msa, truth_msa, col_annotation: str, aa_cmp, msa_annotation_prefix: str = '', scorech=score_chr):
    """
    Compares corresponding columns with each other and writes bool score as column_annotations to test_msa.
    :param test_msa: msa where to save scores
    :param truth_msa: msa against which score will be calculated
    :param col_annotation: name of annotation to save scores
    :param msa_annotation_prefix: if given then writes count and percentage of scores (which tests as True) to annotations of test_sma
    :param aa_cmp: comparator of residues, takes two aa and returns bool
    :return: count and percentage of scores which tests as True
    """

    assert len(test_msa) == len(truth_msa)

    matches = 0
    alnlen = test_msa.get_alignment_length()
    test_msa.column_annotations[col_annotation] = [scorech(False)] * test_msa.get_alignment_length()
    for col in range(alnlen):
        equals = True
        got = None
        for row in range(len(test_msa)):
            tst = test_msa[row].seq[col].upper()
            trh_col = test_msa[row].letter_annotations['fbb_maparray'][col]
            if trh_col == 0:
                equals = False
                break
            
            # substract 1 since maparray starts from 1
            trh = truth_msa[row].seq[trh_col-1].upper()
            if not aa_cmp(tst, trh):
                equals = False
                break
            if got is None:
                got = trh_col
            if got != trh_col:
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


def pairwise(test_msa, truth_msa, col_annotation: str, aa_cmp, msa_annotation_prefix: str= '', scorech=score_chr):
    """
    Compares corresponding cells of corresponding columns with each other and writes scores as column_annotations to test_msa.
    :param test_msa: msa where to save scores
    :param truth_msa: msa against which score will be calculated
    :param col_annotation: name of annotation to save scores
    :param msa_annotation_prefix: if given then writes sum of scores (which tests as True) to annotations of test_sma
    :param aa_cmp: comparator of residues, takes two aa and returns integer score for pair
    :return: sum of scores
    """

    assert len(test_msa) == len(truth_msa)

    sp = 0

    alnlen = test_msa.get_alignment_length()
    test_msa.column_annotations[col_annotation] = [scorech(False)] * test_msa.get_alignment_length()
    for col in range(alnlen):
        colsum = 0
        for row in range(len(test_msa)):
            tst = test_msa[row].seq[col].upper()
            trh_col = test_msa[row].letter_annotations['fbb_maparray'][col]
            if trh_col == 0:
                continue
            
            # substract 1 since maparray starts from 1
            trh = truth_msa[row].seq[trh_col-1].upper()
            score = int(aa_cmp(tst, trh))
            colsum += score
        
        sp += colsum
        test_msa.column_annotations[col_annotation][col] = scorech(colsum)

    test_msa.column_annotations[col_annotation] = ''.join(test_msa.column_annotations[col_annotation])
    if len(msa_annotation_prefix) != 0:
        test_msa.annotations[f'{msa_annotation_prefix}_pw_sum'] = sp
    return sp


def mark_case_by_conservateveness(msa, index_name):
    from Bio.Seq import Seq
    """
    Changes cases of all sequences in msa corresponding to boolean index
    """
    for aln in msa:
        seq = list(aln.seq)
        for i in range(len(seq)):
            if msa.column_annotations[index_name][i] != '1':
                seq[i] = seq[i].lower()
            else:
                seq[i] = seq[i].upper()
        aln.seq = Seq(''.join(seq))
    return msa


def fill_indices(record):
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


def fill_maparray(msa1, msa2, row):
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


def run(args):
    import json
    import pathlib as pl
    from Bio import AlignIO
    from Bio.AlignIO.StockholmIO import StockholmWriter


    def pairwise_scores(x, y):
        if x not in aa_classes or y not in aa_classes:
            return 0
        if x == y:
            return 2
        return 1

    query = AlignIO.read(args.query, args.query_format)
    subject = AlignIO.read(args.subject, args.subject_format)

    for record in query:
        fill_indices(record)
    for record in subject:
        fill_indices(record)
    for i in range(min(len(query), len(subject))):
        fill_maparray(query, subject, i)

    column_wise(query, subject, 'fbb_equals', lambda x, y: x != '-' and x.upper() == y.upper(), msa_annotation_prefix='fbb_equals')
    column_wise(query, subject, 'fbb_similar', are_residues_similar, msa_annotation_prefix='fbb_similar')

    pairwise(query, query, 'fbb_sp', aa_cmp=pairwise_scores, msa_annotation_prefix='fbb_sp')
    thresholder = ConservativenessThresholdMatcher(int(args.threshold), indexes[args.threshold_index])
    if args.mark:
        mark_case_by_conservateveness(query, 'fbb_similar')

    print(json.dumps(query.annotations, indent='  '))
    output = pl.Path(args.output)
    writer = StockholmWriter(output.open('w'))
    writer.pfam_gc_mapping.update({
        'fbb_equals': 'FBB_EQ',
        'fbb_similar': 'FBB_SIM',
        'fbb_sp': 'FBB_SP',
        'fbb_threshold': 'FBB_TH',
    })
    writer.write_alignment(query)


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
        "--mark-conservative-by-case",
        "-m",
        dest='mark',
        help="if given, output sequences will be formated with lower case," +
            "if not conservative column, and upper case otherwise." +
            "uses only pairwise for deciding",
        action='store_true',
    )
    parser.add_argument(
        "--threshold",
        "-d",
        required=True,
        help="threshold for pairwise conservative detection"
    )
    indexes = {
        'eq': 'fbb_equals',
        'sim': 'fbb_similar',
        'sp': 'fbb_sp',
    }
    parser.add_argument(
        "--threshold-index",
        "-i",
        required=True,
        dest='threshold_index',
        help=f"over which index check threshold. supported: {','.join(sorted(indexes.keys()))}"
    )


    args = parser.parse_args()

    if args.threshold_index not in indexes.keys():
        raise ValueError(f'unknown threshold index: {args.threshold_index}')

    run(args)

