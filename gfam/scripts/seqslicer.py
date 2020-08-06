#!/usr/bin/env python
"""Sequence slicer application"""

import sys

from gfam.scripts import CommandLineApp
from gfam.tasks.seqslicer.fasta_fragments_extractor import \
    FastaFragmentsExtractor
from gfam.tasks.seqslicer.slice_file import SliceFile

__authors__ = "Tamas Nepusz, Alfonso E. Romero"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2012, Tamas Nepusz, Alfonso E. Romero"
__license__ = "GPL"


class SeqSlicerApp(CommandLineApp):
    """\
    Usage: %prog [options] [slice_file] [output_file] sequences_file

    Given a sequence database in FASTA format and a list of regions
    to be sliced from those sequences, generates another FASTA file
    that contains the sliced regions.

    slice_file must be a file defining which slices we need. The
    file must be in the following format:

      seqID1 start1 end1
      seqID2 start2 end2
      ...

    When the end position is omitted, the default is the length of
    the sequence. When the start position is also omitted, the whole
    sequence will be returned. Sequence positions are defined by
    integers starting from 1, that is, the first amino acid at the
    N-terminus is at position 1. Both the start and end positions are
    inclusive. If you specify a negative number, this will be interpreted
    as a position from the C-terminus, that is, position -60 is position
    60 counted from the C-terminus.

    sequences_file must be a sequence database in FASTA format.
    """

    short_name = "seqslicer"

    def __init__(self, *args, **kwds):
        super(SeqSlicerApp, self).__init__(*args, **kwds)
        self.output_file = ""

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(SeqSlicerApp, self).create_parser()

        parser.add_option("-a", "--try-alternative-splicing",
                          dest="try_alternative_splicing",
                          default=False,
                          action="store_true",
                          help="try sequenceID.1 if sequenceID is not found")
        parser.add_option("-i", "--ignore-unknown", dest="ignore_unknown",
                          default=False, action="store_true",
                          help="ignore unknown sequence IDs")
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          help="remap sequence IDs using REGEXP",
                          config_key="sequence_id_regexp",
                          dest="sequence_id_regexp")
        parser.add_option("-k", "--keep-ids", dest="keep_ids",
                          action="store_true",
                          help="keep original sequence IDs even"
                               "if this will duplicate "
                               "existing IDs in the output file.")
        return parser

    def run_real(self):
        """Runs the application and returns the exit code"""
        self.output_file = None

        if len(self.args) == 1:
            slice_file = "-"
            seq_file = self.args[0]
        elif len(self.args) == 2:
            slice_file = self.args[0]
            seq_file = self.args[1]
        elif len(self.args) == 3:
            slice_file = self.args[0]
            seq_file = self.args[2]
            self.output_file = self.args[1]
        else:
            self.parser.print_help()
            return 1

        self.log.info(f"Loading slices from {slice_file}...")
        slicer = SliceFile()
        dict_slices = slicer.get_slice_file_as_dict(slice_file)
        self.log.info("Processing fasta file %s...", seq_file)
        extractor = FastaFragmentsExtractor(
            self.output_file,
            self.options.sequence_id_regexp,
            self.options.try_alternative_splicing,
            self.options.keep_ids, self.log)
        return extractor.process_sequences_file(seq_file, dict_slices)


if __name__ == "__main__":
    sys.exit(SeqSlicerApp().run())
