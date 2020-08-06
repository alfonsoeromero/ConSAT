"""This module contains routines for parsing and writing FASTA files.

It is not meant to be a full-fledged FASTA parser (check BioPython_ if you
need one), but it works well in most cases, at least for neatly formatted
FASTA files.
"""

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["Parser", "regexp_remapper", "Writer"]

import re
from itertools import tee
from textwrap import TextWrapper

from gfam.sequence import SeqRecord, Sequence


class Parser(object):
    """Parser for FASTA files.

    Usage example::

        parser = Parser(open("test.ffa"))
        for sequence in parser:
            print sequence.seq

    It also works with remote FASTA files if you use the :mod:`urllib2`
    module::

        from urllib2 import urlopen
        parser = Parser(urlopen("ftp://whatever.org/remote_fasta_file.ffa"))
        for sequence in parser:
            print sequence.seq
    """

    def __init__(self, handle):
        self.handle = handle

    def _lines(self):
        """Iterator that iterates over the interesting lines in a
        FASTA file. This will skip empty lines and trims the remaining
        ones from all unnecessary whitespace. It will also skip the
        lines before the first record enty."""
        handle = self.handle
        while True:
            line = handle.readline()
            if not line:
                return
            if line[0] == ">":
                break

        yield line
        while True:
            line = handle.readline()
            if not line:
                return
            line = line.rstrip().replace("\r", "")
            if not line:
                continue
            yield line

    def sequences(self):
        """Returns a generator that iterates over all the sequences
        in the FASTA file. The generator will yield `SeqRecord`
        objects.
        """
        seq_record, seq = None, []
        for line in self._lines():
            if line[0] == ">":
                # Starting a new sequence here
                if seq_record:
                    seq_record.seq = Sequence("".join(seq))
                    yield seq_record
                descr = line[1:]
                seq_id = descr.split()[0]
                seq_record = SeqRecord(Sequence(), id=seq_id, name=seq_id,
                                       description=descr)
                seq = []
            else:
                # Appending to an existing sequence
                seq.extend(list(line))

        # Don't forget the last sequence
        if seq_record:
            seq_record.seq = Sequence("".join(seq))
            yield seq_record

    def __iter__(self):
        return self.sequences()


def regexp_remapper(iterable, regexp=None, replacement=r'\g<id>'):
    """Regexp-based sequence ID remapper.

    This class takes a sequence iterator as an input and remaps
    each ID in the sequence using a call to `re.sub`.
    `iterable` must yield `SeqRecord` objects; `regexp`
    is the regular expression matching the part to be
    replaced, `replacement` is the replacement string that
    may contain backreferences to `regexp`.

    If `regexp` is ``None`` or an empty string, returns
    the original iterable.
    """
    if not regexp:
        for seq in iterable:
            yield seq
        return

    regexp = re.compile(regexp)
    for seq in iterable:
        seq.id = regexp.sub(replacement, seq.id)
        yield seq


class Writer(object):
    """Writes `SeqRecord` objects in FASTA format to a given file handle."""

    def __init__(self, handle):
        self.handle = handle
        self.wrapper = TextWrapper(width=70)

    def write(self, seq_record):
        """Writes the given sequence record to the file handle passed
        at construction time.
        """
        if seq_record.description:
            print(">%s" % seq_record.description, file=self.handle)
        else:
            print(">%s" % (seq_record.id, ), file=self.handle)
        print("\n".join(self.wrapper.wrap(seq_record.seq)), file=self.handle)


class FastWriter(object):
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
