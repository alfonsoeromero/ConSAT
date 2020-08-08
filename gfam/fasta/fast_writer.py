from itertools import tee


class FastWriter:
    """Writes `SeqRecord` objects in FASTA format to a given file handle.
    This is a version of `Writer` which does not use `TextWrapper` for
    efficiency reasons`"""

    def __init__(self, handle):
        self.handle = handle

    def _pairwise(self, iterable):
        obj1, obj2 = tee(iterable)
        next(obj2, None)
        return zip(obj1, obj2)

    def write(self, seq_record):
        """Writes the given sequence record to the file handle passed
        at construction time
        """
        print(">%s" % seq_record.id, file=self.handle)
        num_residues = len(seq_record.seq)
        pos = list(range(0, num_residues, 70)) + [num_residues]
        for beg_pos, end_pos in self._pairwise(pos):
            print(seq_record.seq[beg_pos:end_pos], file=self.handle)
