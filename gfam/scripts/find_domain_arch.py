#!/usr/bin/env python
"""Application that calculates the domain architecture of each gene and
outputs them in a simple text-based format.
"""
from __future__ import print_function

import sys

from gfam.scripts import CommandLineApp
from gfam.tasks.find_domain_arch.architecture_mapping.architecture_mapping_writer import \
    ArchitectureMappingWriter
from gfam.tasks.find_domain_arch.clustering_file import ClusteringFile

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"


class FindDomainArchitectureApp(CommandLineApp):
    """\
    Usage: %prog [options] interpro_file clustering_file

    Application that calculates the domain architecture of each gene and
    oputs them in a simple text-based format, given the filtered InterPro
    assignments and a clustering of the unknown regions in a separate file.
    """

    short_name = "find_domain_arch"

    def __init__(self, *args, **kwds):
        super(FindDomainArchitectureApp, self).__init__(*args, **kwds)

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(FindDomainArchitectureApp, self).create_parser()
        parser.add_option("-s", "--min-size", dest="min_size",
                          metavar="N",
                          help="consider only clusters with at least N "
                               "different sequences (not just fragments) "
                               "as novel domains",
                          config_key="analysis:find_domain_arch/"
                                     "min_novel_domain_size",
                          default=2, type=int)
        parser.add_option("-S", "--sequences",
                          dest="sequences_file", metavar="FILE",
                          help="FASTA file containing all the sequences "
                               "of the representative gene model",
                          config_key="file.input.sequences", default=None)
        parser.add_option("-i", "--interpro-parent-child-file",
                          dest="interpro_parent_child_file",
                          metavar="FILE",
                          help="use the InterPro parent-child FILE"
                               " to remap IDs",
                          config_key="file.mapping.interpro_parent_child",
                          default=None)
        parser.add_option("--details",
                          dest="details", metavar="FILE",
                          help="print more details about the domain"
                               " architecture into FILE",
                          config_key="generated/"
                                     "file.domain_architecture_details",
                          default=None)
        parser.add_option("--stats",
                          dest="stats", metavar="FILE",
                          help="print genome-level statistics about the domain"
                               " architectures into FILE",
                          config_key="generated/"
                                     "file.domain_architecture_stats",
                          default=None)
        parser.add_option("--new_domains_table",
                          dest="new_domains_table", metavar="FILE",
                          help="prints a table with the new domains,"
                               " one per line, into FILE",
                          config_key="generated/file.new_domains_table",
                          default=None)
        parser.add_option("-n", "--names",
                          dest="interpro_names_file",
                          metavar="FILE",
                          help="use the given FILE to assign InterPro"
                               " IDs to names",
                          config_key="file.mapping.interpro2name",
                          default=None)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          help="remap sequence IDs using REGEXP",
                          config_key="sequence_id_regexp",
                          dest="sequence_id_regexp")
        parser.add_option("--max-overlap", metavar="SIZE",
                          help="sets the maximum overlap size allowed between "
                               "assignments of the same data source."
                               " Default: %default",
                          config_key="max_overlap",
                          dest="max_overlap", type=int, default=20)
        parser.add_option("--previous-table", metavar="FILE", dest="old_table",
                          help="reads a previously built table from a file,"
                               " trying to "
                          "keep the same cluster names as much as possible",
                          config_key="analysis:find_domain_arch/"
                                     "previous_domain_table",
                          default="")
        parser.add_option("-p", metavar="STRING", dest="prefix",
                          help="prefix for the new discovered domains ('NOVEL'"
                               " by default)",
                          config_key="analysis:find_domain_arch/prefix",
                          default='NOVEL')
        return parser

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 2:
            self.error("exactly two input files are expected")

        # create cluster and architecture tables
        clustering_file = ClusteringFile(self.options.min_size,
                                         self.options.old_table,
                                         self.options.prefix)
        mapping = ArchitectureMapping()

        # create app (this prints the output while receiving files
        # and orchestrates everything)

        # print new domain table (from cluster file)
        # print detail file (from architectureassignmenttable)
        # print stats file (from architectureassignmenttable)


if __name__ == "__main__":
    sys.exit(FindDomainArchitectureApp().run())
