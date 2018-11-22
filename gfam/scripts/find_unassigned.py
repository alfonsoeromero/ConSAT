"""Command line script that finds all the regions of a given set
of sequences that are not assigned to any InterPro domain in a
given InterPro domain assignment file"""

from __future__ import print_function
import sys
import shelve
import os
import tempfile

from collections import defaultdict
from gfam import fasta
from gfam.interpro import AssignmentReader
from gfam.scripts import CommandLineApp
from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments
from gfam.utils import open_anything

__authors__ = "Tamas Nepusz, Alfonso E. Romero"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class FindUnassignedApp(CommandLineApp):
    """\
    Usage: %prog [options]

    Given the InterPro assignments for a genome, finds all the regions
    that were not assigned to any InterPro domain and outputs these
    regions in the following format:

      seqID1 start1 end1
      seqID2 start2 end2
      ...

    A sequence ID may appear many times in the first column when multiple
    unassigned regions are present. Starting and ending coordinates are
    both inclusive and start from 1.
    """

    short_name = "find_unassigned"

    def __init__(self, *args, **kwds):
        super(FindUnassignedApp, self).__init__(*args, **kwds)
        self.seqcat = {}
        self.seq_ids_to_length = None
        self.regions = None
        self.low_complexity_regions = None
        self.filename_shelve = ""
        self.sequence_id_regexp = ""

    def create_parser(self):
        """Creates the command line parser used by this script"""
        parser = super(FindUnassignedApp, self).create_parser()
        parser.add_option("-l", "--min-length", dest="min_length",
                          metavar="LENGTH",
                          help="minimum sequence LENGTH needed for a sequence"
                               " in order to include its fragments in "
                               "the output",
                          config_key="analysis:find_unassigned/min_seq_length",
                          default=0, type=int)
        parser.add_option("-f", "--min-fragment-length",
                          dest="min_fragment_length",
                          metavar="LENGTH",
                          help="minimum fragment LENGTH needed in the output",
                          config_key="analysis:find_unassigned/"
                                     "min_fragment_length",
                          default=0, type=int)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          help="remap sequence IDs using REGEXP",
                          config_key="sequence_id_regexp",
                          dest="sequence_id_regexp")
        parser.add_option("-S", "--sequences",
                          dest="sequences_file", metavar="FILE",
                          help="FASTA file containing all the sequences "
                               "of the representative gene model",
                          config_key="analysis:find_unassigned/sequences_file",
                          default=None)
        parser.add_option("--max-overlap", metavar="SIZE",
                          help="sets the maximum overlap size allowed between "
                               "assignments of the same data source."
                               " Default: %default",
                          config_key="max_overlap",
                          dest="max_overlap", type=int, default=20)
        parser.add_option("--low-complexity-regions-file", metavar="FILE",
                          help="file with low complexity regions"
                               "(in segmask format) which will be"
                               " not considered as unassigned",
                          config_key="analysis:find_unassigned/"
                                     "low_complexity_regions_file",
                          dest="low_complexity_file", default=None)
        return parser

    def run_real(self):
        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        if self.options.min_fragment_length < 1:
            self.log.warning("minimum fragment length is not "
                             "positive, assuming 1")
            self.options.min_fragment_length = 1
        self.set_sequence_id_regexp(self.options.sequence_id_regexp)
        self.process_sequences_file_old(self.options.sequences_file)

        self.log.info("Processing assignments file...")
        for infile in (self.args or ["-"]):
            self.process_infile(infile)

        self.log.info("Getting unassigned pieces")
        self.get_unassigned()

        # if there is a file with low complexity regions, process it
        # and remove low complexity regions from unassigned fragments
        if self.options.low_complexity_file:
            self.log.info("Removing low complexity regions")
            self.read_low_complexity_regions(self.options.low_complexity_file)
            self.remove_low_complexity_regions()

        self.log.info("Getting unassigned pieces")
        self.print_unassigned()

    def process_sequences_file(self, fname):
        """ In this version we use `shelve` to save
            memory (the pairs (protein accession, length) are
            stored in a temporary database. See `process_sequences_file_old`
            for the old version.
        """
        self.log.info("Loading sequences from {}...".format(fname))
        parser = fasta.Parser(open_anything(fname))
        parser = fasta.regexp_remapper(parser,
                                       self.sequence_id_regexp)
        self.filename_shelve = os.path.join(tempfile.gettempdir(),
                                            "shelve_file")
        self.seq_ids_to_length = shelve.open(self.filename_shelve)

        for i, seq in enumerate(parser):
            self.seq_ids_to_length[seq.id] = len(seq.seq)
            if i % 1000000 == 0:
                self.log.info("Read {} seqs".format(i))
                self.seq_ids_to_length.sync()
        self.log.info("...loaded")

    def process_sequences_file_old(self, fname):
        """ This is the old version, all the entries are
            loaded into memory
        """
        self.log.info("Loading sequences from %s..." % fname)
        parser = fasta.Parser(open_anything(fname))
        parser = fasta.regexp_remapper(parser,
                                       self.sequence_id_regexp)
        seqs, lens = [], []
        for i, seq in enumerate(parser):
            seqs.append(seq.id)
            lens.append(len(seq.seq))
            if i % 1000000 == 0:
                self.log.info("Read {} seqs".format(i))
        self.log.info("...loaded")
        self.seq_ids_to_length = dict(zip(seqs, lens))

    def process_infile(self, fname, interpro=None):
        self.log.info("Processing input file: %s" % fname)
        for assignment in AssignmentReader(fname):
            try:
                seq = self.seqcat[assignment.id]
            except KeyError:
                seq = SequenceWithAssignments(assignment.id, assignment.length)
                self.seqcat[assignment.id] = seq
            if seq.length != assignment.length:
                raise ValueError("different lengths encountered "
                                 "for %s: %d and %d" % (seq.name,
                                                        seq.length,
                                                        assignment.length))
            if interpro is not None:
                assignment = assignment.resolve_interpro_ids(interpro)
        seq.assign(assignment, interpro= interpro)

    def get_unassigned(self):
        self.regions = []
        for seq_id, seq in self.seqcat.items():
            if seq.length < self.options.min_length:
                continue
            for start, end in seq.unassigned_regions():
                if end-start+1 < self.options.min_fragment_length:
                    continue
                self.regions.append((seq_id, start, end))
        maximum = max(self.options.min_length,
                      self.options.min_fragment_length)
        seqcat_keys = set(self.seqcat.keys())
        for seq_id in set(self.seq_ids_to_length.keys()) - seqcat_keys:
            if self.seq_ids_to_length[seq_id] >= maximum:
                self.regions.append((seq_id, 1,
                                     self.seq_ids_to_length[seq_id]))

    def read_low_complexity_regions(self, file_name):
        self.low_complexity_regions = defaultdict(list)
        current_prot_id = ""

        if self.sequence_id_regexp:
            import re
            regexp = re.compile(self.sequence_id_regexp)
        else:
            regexp = None

        for line in open_anything(file_name):
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                current_prot_id = line.split()[0][1:]
                if regexp is not None:
                    current_prot_id = regexp.sub(r'\g<id>', current_prot_id)
            else:
                (left, _, right) = line.split()
                left = int(left)
                right = int(right)
                self.low_complexity_regions[current_prot_id].append((left,
                                                                     right))

    def substract_interval_from_interval(self, interval1, interval2):
        """Substract from an interval (left, right) another interval
        (leftr,rightr), returning a list of resulting intervals
        (which may be empty) equivalent to the set difference
        """
        (left, right) = interval1
        (leftr, rightr) = interval2

        if leftr > right or rightr < left:
            return [(left, right)]
        elif leftr <= left:
            if rightr >= right:
                return []
            else:
                return [(rightr+1, right)]
        elif rightr >= right:
            # leftr is > left
            return [(left, leftr-1)]
        else:
            # case leftr > left and rightr < right
            return [(left, leftr-1), (rightr+1, right)]

    def substract(self, region, list_to_substract):
        (reg_id, left, right) = region

        # we obtain those regions which overlaps with our region
        # if our region is (left,right) and the set of iterated regions is
        # (leftr, rightr), the overlapping are those (leftr,rightr) where
        subtrahend = [(leftr, rightr) for (leftr, rightr)
                      in list_to_substract
                      if leftr <= right and rightr >= left]

        if not subtrahend:
            # no overlap, the region is returned "as is"
            return [region]
        minuend = [(left, right)]
        for sub in subtrahend:
            minuend.extend(
                self.substract_interval_from_interval(minuend.pop(), sub))
            if not minuend:
                return []
        return [(reg_id, le, ri) for (le, ri) in minuend]

    def remove_low_complexity_regions(self):
        new_regions = []
        for region in self.regions:
            reg_id = region[0]
            if reg_id not in self.low_complexity_regions:
                new_regions.append(region)
            else:
                lcr_ids = self.low_complexity_regions[reg_id]
                new_regions.extend(self.substract(region, lcr_ids))
        self.regions = new_regions

    def print_unassigned(self):
        maximum = max(self.options.min_length,
                      self.options.min_fragment_length)
        for region in self.regions:
            if region[2] - region[1] + 1 >= maximum:
                print("%s\t%d\t%d" % region)

    def set_sequence_id_regexp(self, regexp):
        self.sequence_id_regexp = regexp


if __name__ == "__main__":
    sys.exit(FindUnassignedApp().run())
