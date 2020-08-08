from gfam.sequence import SeqRecord, Sequence


class Parser:
    """Parser for FASTA files.

    Usage example::

        parser = Parser(open("test.ffa"))
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
