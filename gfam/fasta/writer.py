from textwrap import TextWrapper


class Writer:
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
