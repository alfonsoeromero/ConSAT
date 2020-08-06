"""Command line script that finds all the regions of a given set
of sequences that are not assigned to any InterPro domain in a
given InterPro domain assignment file"""

import sys

from gfam.scripts import CommandLineApp
from gfam.tasks.find_unassigned.task import FindUnassignedTask

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

    def _filter_parameters(self) -> None:
        if self.options.min_fragment_length < 1:
            self.log.warning("minimum fragment length is not "
                             "positive, assuming 1")
            self.options.min_fragment_length = 1

    def run_real(self):
        self._filter_parameters()
        assignment_file = (self.args or ["-"])[0]

        task = FindUnassignedTask(self.options.max_overlap,
                                  self.options.min_fragment_length,
                                  self.options.min_length,
                                  self.options.sequence_id_regexp,
                                  self.log)

        task.run(assignment_file,
                 self.options.sequences_file,
                 self.options.low_complexity_file)


if __name__ == "__main__":
    sys.exit(FindUnassignedApp().run())
