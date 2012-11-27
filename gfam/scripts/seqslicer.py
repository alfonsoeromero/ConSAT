#!/usr/bin/env python
"""Sequence slicer application"""

import sys

from gfam import fasta
from gfam.sequence import SeqRecord
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything

__authors__  = "Tamas Nepusz, Alfonso E. Romero"
__email__   = "tamas@cs.rhul.ac.uk"
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
        self.seqs = None

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(SeqSlicerApp, self).create_parser()

        parser.add_option("-a", "--try-alternative-splicing",
                dest="try_alternative_splicing", default=False,
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
                help="keep original sequence IDs even if this will duplicate "
                     "existing IDs in the output file.")
        return parser

    def load_slice_file(self, slice_file):
        """Loads the slice file into a dictionary of lists"""
        self.log.info("Loading slices from %s..." % slice_file)

        self.parts = defaultdict()

        for line in open_anything(slice_file):
            parts = line.strip().split()
            if not parts:
                continue
            seq_id = parts[0]
            (left, right) = (1, None)

            if len(parts) == 3:
                # Three cases: (a) both limits (left,right) are specified
                (left, right) = (int(parts[1]), int(parts[2]))
            elif len(parts) == 2:
                # (b) only the left limit is specified
                left = int(parts[1])
                # (c) neither left nor right limits are specified 
                # (this is why we se left=1 before)

            self.parts[seq_id].append((left, right))

    def process_sequences_file(self, seq_file):
        """Processes the sequences one by one, extracting all the pieces into
        an output fasta file"""
        self.log.info("Processing fasta file %s..." %seq_file)

        parser = fasta.Parser(open_anything(seq_file))
        parser = fasta.regexp_remapper(parser, 
            self.options.sequence_id_regexp)

        ids_to_process = set(self.parts.keys())

        writer = fasta.Writer(sys.stdout)
        if self.output_file is not None:
            self.output_fd = open(output_file,"r")
            writer_file = fasta.Writer(output_fd)

        for seq in parser:
            seq_id = seq.id
            if seq_id not in self.parts:
                if self.options.try_alternative_splicing:
                    seq_id = seq_id.strip().rstrip(".1")
                    if seq_id not in self.parts:
                        continue
                else:
                    continue
 
            sequence = seq.seq
            length_seq = len(sequence)
            ids_to_process.remove(seq_id)

            for (left, right) in self.parts[seq_id]:

                if left == 0:
                    self.log.warning("Ignoring fragment ID: %s, "
                        "requested start position is zero" % seq_id)
                    continue
                elif left < 0:
                    left = length_seq + left + 1

                if right is None:
                    right = length_seq
                elif right < 0:
                    right = length_seq + right + 1
                elif right > length_seq:
                    #just in case...
                    right = length_seq
                elif right == 0:
                    self.log.warning("Ignoring fragment ID: %s, "
                        "requested end position is zero" %seq_id)
                    continue

                if left >= right:
                    #again, just in case
                    self.log.warning("Problem with fragment of %s, "
                        "the right part is smaller than the left" % seq_id)
                    continue

                if not self.options.keep_ids:
                    new_id = "%s:%d-%d" % (seq_id, left, right)
                else:
                    new_id = seq_id

                new_record = SeqRecord(sequence[(left-1):right],
                        id=new_id, name=seq.name, description="")
                writer.write(new_record)

	        if self.output_file is not None:
                    writer_file.write(new_record)

        if self.output_file is not None:
            output_fd.close()

        if len(ids_to_process) > 0:
            self.log.fatal("The following identifiers of sequences (%s) were"
                    "found in the fragments file, but not in the fasta file"
                    % ",".join(ids_to_process))
            return 1

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

        self.load_slice_file(slice_file)
        self.process_sequences_file(seq_file)

if __name__ == "__main__":
    sys.exit(SeqSlicerApp().run())

