#!/usr/bin/env python

import sys

from gfam.scripts import CommandLineApp
from gfam.tasks.extract_gene_ids.task import ExtractGeneIDsTask

__author__ = "Tamas Nepusz, Alfonso E. Romero"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2020, Tamas Nepusz, Alfonso E. Romero"
__license__ = "GPL"


class ExtractGeneIDsApp(CommandLineApp):
    """\
    Usage: %prog [options] [result_file]

    Extracts the gene IDs from a FASTA file.
    """

    def __init__(self):
        self.task = ExtractGeneIDsTask()

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(ExtractGeneIDsApp, self).create_parser()
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          help="remap sequence IDs using REGEXP",
                          config_key="sequence_id_regexp",
                          dest="sequence_id_regexp")
        return parser

    def run_real(self):
        """Runs the application"""
        if not self.args:
            infiles = ["-"]
        else:
            infiles = self.args

        for infile in infiles:
            self.task.process_file(infile,
                                   self.options.sequence_id_regexp)


if __name__ == "__main__":
    sys.exit(ExtractGeneIDsApp().run())
