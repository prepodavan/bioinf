from collections.abc import Iterable


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


def residues_from_class(classes: Iterable[str]) -> Iterable[str]:
    cls = [cls.upper() for cls in classes]
    return set([it[0] for it in filter(lambda it: it[1] in cls, aa_classes.items())])


class Column:
    def __init__(self, pos: int, rows: int) -> None:
        self.pos = pos
        self.rows = rows
        self.residues = dict()
        self.classes  = dict()

    def __str__(self) -> str:
        indices  = sum(self.residues.values())
        residues = ''.join(self.residues.keys())
        s = f'{self.pos + 1}:[{residues}]'
        if indices == self.rows:
            return s
        return f'{s}[{indices}/{self.rows}]'

    def add(self, residue: str) -> None:
        residue = residue.upper()
        self.residues[residue] = 1 + self.residues.get(residue, 0)
        cls = aa_classes.get(residue)
        if cls is not None:
            self.classes[cls] = 1 + self.classes.get(cls, 0)

    def remove(self, residues: Iterable[str]) -> int:
        """
        Removes all residues from column
        :param residues: set of residues to delete
        :return: number of deleted residues
        """
        deleted = 0
        for res in residues:
            if res not in self.residues.keys():
                continue
            deleted += 1
            num = self.residues.pop(res)
            cls = aa_classes.get(res)
            if cls is None or cls not in self.classes.keys():
                continue
            self.classes[cls] -= num
            if self.classes[cls] == 0:
                self.classes.pop(cls)
        return deleted


def cons_to_md(alnrows: int, full: Iterable[Column], func: Iterable[Column], reliable: Iterable[Iterable[int]]) -> str:
    import io


    maxlen = max([len(str(col)) for col in full + func + ['seed', '100%-cons', 'func-cons', 'max-reliable']])
    header_cell = '-' * maxlen
    buf = io.StringIO()
    buf.write(f"|{'seed'.ljust(maxlen)}|{'100%-cons'.ljust(maxlen)}|{'func-cons'.ljust(maxlen)}|{'max-reliable'.ljust(maxlen)}|\n" + \
        "|" + '|'.join([header_cell]*4) + "|\n",
    )
    rows = max(len(full), len(func), len(reliable))
    for i in range(rows):
        buf.write('|' + f'{alnrows if i == 0 else ""}'.ljust(maxlen) + '|')
        for var in [full, func, reliable]:
            buf.write(f'{var[i] if i < len(var) else ""}'.ljust(maxlen) + '|')
        buf.write(f'\n')

    return buf.getvalue()


def block_contains(
        msa,
        col_start: int,
        col_end: int,
        symbol: str = '-',
        row_start: int = 0,
        row_end: int = -1,
) -> bool:
    """
    Check whether block specified by coordinates contains symbol
    :param row_end: if not given or given negative, then uses len(msa)
    """
    if row_end < 0:
        row_end = len(msa)
    assert 0 <= col_start < col_end and 0 <= row_start < row_end <= len(msa)

    for row in range(row_start, row_end):
        seq = msa[row].seq
        assert col_end <= len(seq)
        if symbol in seq[col_start:col_end]:
            return True

    return False


def find_conservatives(
        msa,
        full_threshold: int,
        func_threshold: int,
        mark: str = None,
    ):
    full = list()
    func = list()
    reliable = list()

    alnlen = msa.get_alignment_length()
    msa.column_annotations['FBB_FULL'] = ['-'] * alnlen
    msa.column_annotations['FBB_FUNC'] = ['-'] * alnlen
    msa.column_annotations['FBB_RLBL'] = ['-'] * alnlen

    for col in range(len(msa[0].seq)):
        column = Column(col, len(msa))

        for row in range(len(msa)):
            assert len(msa[row].seq) == len(msa[0].seq)
            column.add(msa[row].seq[col])

        if len(column.classes) == 0 or len(column.residues) == 0:
            continue

        if '-' in column.residues.keys():
            continue

        if len(column.residues) == 1:
            full.append(column)
            msa.column_annotations['FBB_FULL'][col] = '1'
            continue

        full_cons = list(filter(lambda x: x[1] >= full_threshold, column.residues.items()))
        full_cons.sort(key=lambda x: -1 * x[1])
        if len(full_cons) != 0:
            column.remove(set(column.residues.keys()) - {full_cons[0][0]})
            full.append(column)
            msa.column_annotations['FBB_FULL'][col] = '1'
            continue

        if len(column.classes) == 1:
            func.append(column)
            msa.column_annotations['FBB_FUNC'][col] = '1'
            continue

        func_cons = list(filter(lambda x: x[1] >= func_threshold, column.classes.items()))
        func_cons.sort(key=lambda x: -1 * x[1])
        if len(func_cons) != 0:
            column.remove(set(column.residues.keys()) - residues_from_class(func_cons[0][0]))
            func.append(column)
            msa.column_annotations['FBB_FUNC'][col] = '1'

    msa.column_annotations['FBB_FULL'] = ''.join(msa.column_annotations['FBB_FULL'])
    msa.column_annotations['FBB_FUNC'] = ''.join(msa.column_annotations['FBB_FUNC'])

    if len(full) + len(func) < 2:
        msa.column_annotations['FBB_RLBL'] = ''.join(msa.column_annotations['FBB_RLBL'])
        if mark is not None and len(mark) != 0:
            mark_case_by_conservateveness(msa, mark)
        return msa, full, func, reliable

    cols = sorted(full + func, key=lambda x: x.pos)
    flag = False
    start = 0
    for i in range(1, len(cols)):
        gaped = block_contains(msa, cols[i-1].pos, cols[i].pos + 1, '-')
        if not flag and not gaped:
            flag = True
            start = cols[i-1].pos
        elif flag and gaped:
            flag = False
            reliable.append((start + 1, cols[i-1].pos + 1))
            for col in range(start, cols[i-1].pos + 1):
                msa.column_annotations['FBB_RLBL'][col] = '1'

    msa.column_annotations['FBB_RLBL'] = ''.join(msa.column_annotations['FBB_RLBL'])
    if mark is not None and len(mark) != 0:
        mark_case_by_conservateveness(msa, mark)
    return msa, full, func, reliable


def mark_case_by_conservateveness(msa, index_name: str):
    """
    Changes cases of all sequences in msa corresponding to boolean index
    """
    from Bio.Seq import Seq


    for aln in msa:
        seq = list(aln.seq)
        for i in range(len(seq)):
            if msa.column_annotations[index_name][i] != '1':
                seq[i] = seq[i].lower()
            else:
                seq[i] = seq[i].upper()
        aln.seq = Seq(''.join(seq))
    return msa


def write_to_sth(msa, output):
    from Bio.AlignIO.StockholmIO import StockholmWriter


    writer = StockholmWriter(output)
    for annotation_name in msa.column_annotations.keys():
        if annotation_name.startswith('FBB_'):
            writer.pfam_gc_mapping[annotation_name] = annotation_name
    writer.write_alignment(msa)


def find_conservatives_from_file(
        input_filepath: str = None,
        input_format: str = None,
        full_threshold: int = None,
        func_threshold: int = None,
        output: str = None,
        mark: str = None,
    ):
    from Bio import AlignIO


    msa, full, func, reliable = find_conservatives(
        AlignIO.read(input_filepath, input_format),
        int(full_threshold),
        int(func_threshold),
        mark,
    )

    if output is None or len(output) == 0:
        return msa, full, func, reliable

    with open(output, 'w') as f:
        write_to_sth(msa, f)

    return msa, full, func, reliable


if __name__ == '__main__' and '__file__' in globals():
    import argparse


    parser = argparse.ArgumentParser(
        description='This program searches for convservative, functionally convservative and maximally reliable columns of MSA.',
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        dest='input_filepath',
        help="input filepath with alignment",
    )
    parser.add_argument(
        "--input-format",
        "-f",
        required=True,
        dest='input_format',
        help="format of query. required by Bio.AlignIO",
    )
    parser.add_argument(
        "--full-threshold",
        "-l",
        required=True,
        dest='full_threshold',
        help="threshold for detecting fully conservative column"
    )
    parser.add_argument(
        "--func-threshold",
        "-n",
        required=True,
        dest='func_threshold',
        help="threshold for detecting fully conservative column"
    )
    parser.add_argument(
        "--output",
        "-o",
        default='',
        help="filepath to write down stockholm output. print markdown table to stdout, if not given",
    )
    parser.add_argument(
        "--mark-conservative-by-case",
        "-m",
        dest='mark',
        default='',
        help="index name which will store sequences formated with lower case," +
            "if not conservative column, and upper case otherwise." +
            "stores nothing if not given",
    )

    args = parser.parse_args()
    msa, full, func, reliable = find_conservatives_from_file(**vars(args))
    if args.output == '':
        print(cons_to_md(len(msa), full, func, reliable))

