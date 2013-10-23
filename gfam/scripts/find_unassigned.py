#!/usr/bin/env python
"""Command line script that finds all the regions of a given set
of sequences that are not assigned to any InterPro domain in a
given InterPro domain assignment file
"""

import bisect
import optparse
import sys

from collections import defaultdict
from gfam import fasta
from gfam.interpro import AssignmentReader
from gfam.scripts import CommandLineApp
from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments
from gfam.utils import open_anything

__authors__  = "Tamas Nepusz, Alfonso E. Romero"
__email__   = "tamas@cs.rhul.ac.uk"
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

    def create_parser(self):
        """Creates the command line parser used by this script"""
        parser = super(FindUnassignedApp, self).create_parser()
        parser.add_option("-l", "--min-length", dest="min_length",
                metavar="LENGTH",
                help="minimum sequence LENGTH needed for a sequence in order to include its fragments in the output",
                config_key="analysis:find_unassigned/min_seq_length",
                default=0, type=int)
        parser.add_option("-f", "--min-fragment-length", dest="min_fragment_length",
                metavar="LENGTH",
                help="minimum fragment LENGTH needed in the output",
                config_key="analysis:find_unassigned/min_fragment_length",
                default=0, type=int)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP",
                config_key="sequence_id_regexp",
                dest="sequence_id_regexp")
        parser.add_option("-S", "--sequences",
                dest="sequences_file", metavar="FILE",
                help="FASTA file containing all the sequences of the representative gene model",
                config_key="analysis:find_unassigned/sequences_file",
                default=None)
        parser.add_option("--max-overlap", metavar="SIZE",
                help="sets the maximum overlap size allowed between "
                     "assignments of the same data source. Default: %default",
                config_key="max_overlap",
                dest="max_overlap", type=int, default=20)
        parser.add_option("--low-complexity-regions-file", metavar="FILE",
                help="file with low complexity regions (in segmask format) which will be not considered as unassigned",
                config_key="analysis:find_unassigned/low_complexity_regions_file",
                dest="low_complexity_file",default=None)
        return parser

    def run_real(self):
        AssignmentOverlapChecker.max_overlap = self.options.max_overlap

        if self.options.min_fragment_length < 1:
            self.log.warning("minimum fragment length is not positive, assuming 1")
            self.options.min_fragment_length = 1
        self.set_sequence_id_regexp(self.options.sequence_id_regexp)
        self.process_sequences_file(self.options.sequences_file)

        for infile in (self.args or ["-"]):
            self.process_infile(infile)

        self.get_unassigned()

        # if there is a file with low complexity regions, process it
        # and remove low complexity regions from unassigned fragments
        if self.options.low_complexity_file:
            self.read_low_complexity_regions(self.options.low_complexity_file)
            self.remove_low_complexity_regions()

        self.print_unassigned()

    def process_sequences_file(self, fname):
        self.log.info("Loading sequences from %s..." % fname)
        self.seq_ids_to_length = {}
        parser = fasta.Parser(open_anything(fname))
        parser = fasta.regexp_remapper(parser,
                self.sequence_id_regexp
        )
        for seq in parser:
            self.seq_ids_to_length[seq.id] = len(seq.seq)

    def process_infile(self, fname, interpro=None):
        self.log.info("Processing input file: %s" % fname)
        import sys
        for assignment in AssignmentReader(fname):
            try:
                seq = self.seqcat[assignment.id]
            except KeyError:
                seq = SequenceWithAssignments(assignment.id, assignment.length)
                self.seqcat[assignment.id] = seq
            if seq.length != assignment.length:
                raise ValueError, "different lengths encountered for %s: %d and %d" % (seq.name, seq.length, assignment.length)
            if interpro is not None:
                assignment = assignment.resolve_interpro_ids(interpro)
            seq.assign(assignment)

    def get_unassigned(self):
        self.regions = []
        for seqID, seq in self.seqcat.iteritems():
            if seq.length < self.options.min_length:
                continue
            for start, end in seq.unassigned_regions():
                if end-start+1 < self.options.min_fragment_length:
                    continue
                self.regions.append((seqID, start, end))
        maximum = max(self.options.min_length, self.options.min_fragment_length)
        for seqID in set(self.seq_ids_to_length.keys()) - set(self.seqcat.keys()):
            if self.seq_ids_to_length[seqID] >= maximum:
                self.regions.append((seqID, 1, self.seq_ids_to_length[seqID]))

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
                self.low_complexity_regions[current_prot_id].append((int(left),int(right)))


    def substract_interval_from_interval(self, interval1, interval2):
        """Substract from an interval (left, right) another interval 
        (leftr,rightr), returning a list of resulting intervals 
        (which may be empty) equivalent to the set difference
        """
        (left,   right) = interval1
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
            return [(left,leftr-1), (rightr+1, right)]

    def substract(self, region, list_to_substract):
        (id, left, right) = region

        # we obtain those regions which overlaps with our region
        # if our region is (left,right) and the set of iterated regions is
        # (leftr, rightr), the overlapping are those (leftr,rightr) where
        # 
        subtrahend = [(leftr, rightr) for (leftr, rightr) \
                in list_to_substract \
                if leftr <= right and rightr >= left]

        if len(subtrahend) == 0:
            # no overlap, the region is returned "as is"
            return [region]
        else :
            minuend = [(left,right)]
            for sub in subtrahend:
                minuend.extend( \
                    self.substract_interval_from_interval(minuend.pop(), sub))
                if len(minuend) == 0:
                    return []
            return [(id, le, ri) for (le, ri) in minuend]

    def remove_low_complexity_regions(self):
        new_regions = []
        for region in self.regions:
            (id, left, right) = region
            if id not in self.low_complexity_regions:
                new_regions.append(region)
            else:
                new_regions.extend(self.substract(region, self.low_complexity_regions[id])) 

        self.regions = new_regions

    def print_unassigned(self):
        maximum = max(self.options.min_length, self.options.min_fragment_length)
        for region in self.regions:
            if region[2] - region[1] + 1 >= maximum:
                print "%s\t%d\t%d" % region

    def set_sequence_id_regexp(self, regexp):
        self.sequence_id_regexp = regexp


if __name__ == "__main__":
    sys.exit(FindUnassignedApp().run())
