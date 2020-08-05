#!/usr/bin/env python

from __future__ import print_function

import sys
from typing import List, Set

from gfam.assignment import AssignmentOverlapChecker, EValueFilter
from gfam.scripts import CommandLineApp
from gfam.tasks.assignment_source_filter.assignment_filter import \
    AssignmentFilter
from gfam.tasks.assignment_source_filter.assignment_reader_with_filters import\
    AssignmentReaderWithFilters
from gfam.tasks.assignment_source_filter.assignment_source_filter_task import \
    AssignmentSourceFilterTask
from gfam.tasks.assignment_source_filter.exclusion_logger.exclusion_logger\
    import ExclusionLogger
from gfam.tasks.assignment_source_filter.interpro_file_factory import \
    InterproFileFactory
from gfam.tasks.assignment_source_filter.stages_from_config_reader import \
    StagesFromConfigReader
from gfam.tasks.assignment_source_filter.valid_sequence_ids_factory import \
    ValidSequenceIdsFactory

__author__ = "Tamas Nepusz, Alfonso E. Romero"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2020, Tamas Nepusz, Alfonso E. Romero"
__license__ = "GPL"


class AssignmentSourceFilterApp(CommandLineApp):
    """\
    Usage: %prog [options] [assignment_file]

    Tries to determine a consensus domain architecture for sequences based
    on the output of IPRScan. Domains that have a corresponding InterPro ID
    have priority over those who don't have one.

    This program expects incoming assignments from the standard input or from
    a given file and print the selected ones to the standard output.

    The assignment process has three stages by default; stages can be
    configured in the config file. The default setup is as follows:

        1. For a given sequence, take the source having the maximal coverage
           and use it as a primary assignment. In this step, HMMPanther and
           Gene3D domains are excluded.

        2. Loop over the unused domains and try to augment the primary
           assignment with them. Domain insertions and overlaps are allowed
           only if both
           domains (the one being inserted and the one which is already
           there in the assignment) have the same data source.
           In this step, HMMPanther and Gene3D domains are still excluded.

        3. Try step 2 again with HMMPanther and Gene3D domains.
    """

    short_name = "assignment_source_filter"

    def create_parser(self):
        """Creates the command line parser used by this script"""
        parser = super(AssignmentSourceFilterApp, self).create_parser()
        parser.add_option("-x", "--exclude", dest="ignored",
                          metavar="SOURCE",
                          help="add SOURCE to the list of ignored sources",
                          config_key="analysis:iprscan_filter/" +
                                     "untrusted_sources",
                          action="append", default=[])
        parser.add_option("-e", "--e-value", dest="max_e",
                          metavar="THRESHOLD",
                          help="E-value THRESHOLD to filter assignments",
                          config_key="analysis:iprscan_filter/" +
                                     "e_value_thresholds",
                          default="inf")
        parser.add_option("-i", "--interpro-file", dest="interpro_file",
                          metavar="FILE",
                          help="use the InterPro parent-child "
                               " FILE to remap IDs",
                          config_key="analysis:iprscan_filter/" +
                                     "interpro_parent_child_mapping",
                          default=None)
        parser.add_option("-g", "--gene-ids", dest="gene_id_file",
                          metavar="FILE", help="only consider those IDs which "
                          "are present in the list in the given FILE",
                          config_key="generated/file.valid_gene_ids",
                          default=None)
        parser.add_option("--log-exclusions", dest="exclusions_log_file",
                          metavar="FILE",
                          help="log excluded sequences to the given FILE "
                               "for debugging purposes",
                          default=None,
                          config_key="DEFAULT/file.log.iprscan_exclusions")
        parser.add_option("--max-overlap", metavar="SIZE",
                          help="sets the maximum overlap size allowed between "
                               "assignments of the same data source. "
                               "Default: %default",
                          config_key="max_overlap",
                          dest="max_overlap", type=int, default=20)
        return parser

    def _read_ignored_sources(self, sources: List[str]) -> Set[str]:
        ignored = set()
        for ignored_source in sources:
            parts = ignored_source.split()
            ignored.update(parts)
        return ignored

    def _check_args(self) -> None:
        if not self.args:
            self.args = ["-"]
        elif len(self.args) > 1:
            self.error("Only one input file may be given")

    def run_real(self):
        """Runs the application"""
        # sets the algorithm main parameters
        AssignmentOverlapChecker.max_overlap = self.options.max_overlap
        interpro = InterproFileFactory.get_from_file(
            self.options.interpro_file)

        exclusion_logger = ExclusionLogger(
            self.options.exclusions_log_file)
        valid_sequence_ids = ValidSequenceIdsFactory(
            self.log, self.options.gene_id_file).get()
        ignored = self._read_ignored_sources(self.options.ignored)

        self._check_args()
        stages_from_config = StagesFromConfigReader(
            self.parser).get_stages_from_config()
        evalue_filter = EValueFilter.from_string(self.options.max_e)
        assignment_reader_with_filters = AssignmentReaderWithFilters(
            self.args[0], evalue_filter, ignored)
        assignment_filter = AssignmentFilter(exclusion_logger,
                                             stages_from_config,
                                             interpro,
                                             self.log)
        task = AssignmentSourceFilterTask(assignment_reader_with_filters,
                                          assignment_filter,
                                          self.log)
        task.process_infile(valid_sequence_ids)
        exclusion_logger.close()


if __name__ == "__main__":
    sys.exit(AssignmentSourceFilterApp().run())
