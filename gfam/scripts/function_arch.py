#!/usr/bin/env python
"""Application that transfers the function of a set of proteins A (from
a GOA file) whose architecture has been computed, from a set of proteins
B, whose architecture is provided.

The GO terms transferred are printed, for each protein id, in a text
readable format
"""

from collections import defaultdict
from gfam.go import Tree as GOTree
from gfam.go.overrepresentation import OverrepresentationAnalyser
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything
from gfam.utils import bidict
from gfam.architecture import ArchitectureFileReaderPerArch as ArchReader
import operator
import optparse
import sys
import os

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"


class TransferFunctionFromDomainArch(CommandLineApp):
    """\
    Usage: %prog [options] gene_ontology_file architecture_file_A goa_B

    Application that transfers functions to a set of proteins A (from
    a GOA file) whose architecture has been computed, from a set of proteins
    B, whose architecture is provided. The transference is performed
    architecture-wise and overrepresentation analysis is carried out on
    the set of assigned labels.
    The final results are printed to the standard output 
    """

    short_name = "function_arch"

    def __init__(self, *args, **kwds):
        super(TransferFunctionFromDomainArch, self).__init__(*args, **kwds)
        self.evidence = {"EXP": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                 "TAS", "IC"],
                         "ALL_BUT_IEA": ["EXP", "IDA", "IPI", "IMP", "IGI",
                                         "IEP", "TAS", "IC", "ISS", "ISO",
                                         "ISA", "ISM", "IGC", "IBA", "IBD",
                                         "IKR", "IRD", "RCA"],
                         "ALL": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                 "TAS", "IC", "ISS", "ISO", "ISA", "ISM",
                                 "IGC", "IBA", "IBD", "IKR", "IRD", "RCA",
                                 "IEA"]}

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(TransferFunctionFromDomainArch, self).create_parser()
        parser.add_option("-m", "--min-coverage", dest="minimum_coverage",
                          type=float, default=80.0, metavar="PERCENTAGE",
                          config_key="analysis:function_arch/minimum_coverage",
                          help="Minimum % of coverage allowed for an " +
                               "architecture to transfer/receive function")
        parser.add_option("-p", "--p-value", dest="max_pvalue",
                          type=float, default=0.05, metavar="FLOAT",
                          config_key="analysis:function_arch/max_pvalue",
                          help="Maximum allowed p-value for " +
                          "overrepresented GO terms")
        parser.add_option("-e", "--ev_codes", dest="ev_codes", 
                          default="EXP", metavar="EXP/ALL_BUT_IEA/ALL",
                          config_key="analysis:function_arch/ev_codes",
                          choices=("EXP", "ALL_BUT_IEA", "ALL"),
                          help="Evidence codes to use: EXP (experimental), " +
                               "ALL_BUT_IEA (all except electronic, IEA," +
                               "NAS and ND), ALL (all evidence codes)")
        parser.add_option("-s", "--source_arch", dest="source_arch",
                          metavar="FILE", config_key="file.function.domain_arch_table",
                          help="Architecture file to transfer from. If ommited the " +
                          "same architecture file will be used on itself")
        parser.add_option("-b", "--both", dest="both",
                         config_key="analysis:function_arch/transfer_from_both",
                         help="transfer both from the given architectures "
                         " and the computed ones (GOA file is assummed to have "
                         " both protein sets")
        parser.add_option("-a", "--arch_file", dest="results_by_arch",
                         metavar="results_by_arch", 
                         config_key="generated/file.function_arch.general_arch_file",
                         help="File where the results per architecture will"
                         "be written (optional)")
        return parser

    def read_goa_file(self, goa_file, ev_codes):
        """Reads the GOA file, return a defaultdict that,
        for each protein id it has a set of strings with the
        associated GO terms
        """
        d = bidict()
        for line in open_anything(goa_file):
            if not line.startswith(("!", "#")):
                # split line, obtain protein_id, go_term, and ev_code
                fields = line.split("\t", 7)
                prot_id, goterm, evcode = fields[1], fields[4], fields[6]

                if evcode in ev_codes:
                    d.add_left(prot_id, self.go_tree.lookup(goterm))
        self.log.info("GOA file read. " + str(d.len_left()) + " proteins loaded")
        return d

    def _transfer_from_same_file(self, goa, arch_file):
        """ Transfer function from same architecture file
        """
        ora = OverrepresentationAnalyser(self.go_tree, goa,
                                         confidence=self.options.max_pvalue,
                                         min_count=1, correction='None')
        cov = self.options.minimum_coverage / 100.0
        self.log.info("Transferring function from same file. Min coverage=" + str(cov))
        all_annotated = frozenset(goa.left.keys())
        if self.options.results_by_arch:
            out = open(self.options.results_by_arch, "w")
        for arch, prots in ArchReader(arch_file, cov):
            targets = set(prots)
            annotated_prots = targets & all_annotated
            if not annotated_prots or arch == "NO_ASSIGNMENT":
                # if there is no annotation for proteins in the arch...
                print (os.linesep * 2).join(prots), os.linesep
                if self.options.results_by_arch:
                    out.write(arch + "\n\n")
                continue

            lines = os.linesep.join(["  %.4f: %s (%s)" % 
                (p_value, term.id, term.name)
                for term, p_value in ora.test_group(targets)])
            if self.options.results_by_arch:
                out.write(arch + "\n")
                out.write(lines + "\n")
                out.write("\n")
            for prot in prots:
                print prot
                if prot in annotated_prots:
                    for term, p_value in ora.test_group(targets - set([prot])):
                        print "  %.4f: %s (%s)" % (p_value, term.id, term.name)
                else:
                    print lines
                print
        if self.options.results_by_arch:
            out.close()

    def _transfer_from_both(self, goa, arch_target, arch_source):
        """ Transfer from both an external file and the same file using a single GOA 
        file
        """
        ora = OverrepresentationAnalyser(self.go_tree, goa,
                                         confidence=self.options.max_pvalue,
                                         min_count=1, correction='None')
        cov = self.options.minimum_coverage / 100.0
        self.log.info("Transferring function from both files. Min coverage=" + str(cov))
        all_annotated = frozenset(goa.left.keys()) # all annotated proteins
        prots_per_arch = dict()
        if self.options.results_by_arch:
            out = open(self.options.results_by_arch, "w")

        self.log.info("\t Source architecture: " + arch_source)
        self.log.info("\t Target (and source, as well) architecture: " + arch_target)

        for arch, prots in ArchReader(arch_source, cov):
            prots_per_arch[arch] = prots
        for arch, prots in ArchReader(arch_target, cov):
            other_prots = set()
            if arch in prots_per_arch:
                other_prots = prots_per_arch[arch]
            annotated_prots = (other_prots | prots) & all_annotated
            if not annotated_prots or arch == "NO_ASSIGNMENT":
                # if there is no annotation for proteins in the arch...
                print (os.linesep * 2).join(prots), os.linesep
                if self.options.results_by_arch:
                    out.write(arch + "\n\n")
                continue
            lines = os.linesep.join(["  %.4f: %s (%s)" % 
                (p_value, term.id, term.name)
                for term, p_value in ora.test_group(targets)])
            if self.options.results_by_arch:
                out.write(arch + "\n")
                out.write(lines)
                out.write("\n")
            for prot in prots:
                print prot
                if prot in annotated_prots:
                    for term, p_value in ora.test_group(targets - set([prot])):
                        print "  %.4f: %s (%s)" % (p_value, term.id, term.name)
                else:
                    print lines
                print
        if self.options.results_by_arch:
            out.close()


    def _transfer_from_other_file(self, goa, arch_target, arch_source):
        ora = OverrepresentationAnalyser(self.go_tree, goa,
                                         confidence=self.options.max_pvalue,
                                         min_count=1, correction='None')
        cov = self.options.minimum_coverage / 100.0
        goterms = dict()
        for arch, prots in ArchReader(arch_source, cov):
            if arch != "NO_ASSIGNMENT":
                goterms[arch] = ora.test_group(prots)

        if self.options.results_by_arch:
            with open(self.options.results_by_arch, "w") as out:
                for arch in sorted(goterms.keys()):
                    out.write(arch + "\n")
                    for term, p_value in goterms[arch]:
                        line = "  %.4f: %s (%s)" % (p_value, term.id, term.name)
                        out.write(line)
                        out.write("\n")
                    out.write("\n")

        for arch, prots in ArchReader(arch_target, cov):
            if arch in goterms:
                for prot in prots:
                    print prot
                    for term, p_value in goterms[arch]:
                        print "  %.4f: %s (%s)" % (p_value, term.id, term.name)
                    print
            else:
                for prot in prots:
                    print prot
                    print

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 3:
            self.error("exactly three input files are expected" + str(len(self.args)))

        if self.options.ev_codes not in self.evidence:
            self.error("The three valid types of evidence codes are: " +
                       "EXP, ALL_BUT_IEA, and ALL")

        if self.options.max_pvalue < 0.0 or self.options.max_pvalue > 1.0:
            self.error("The maximum p-value should be between 0.0 and 1.0")

        codes = set(self.evidence[self.options.ev_codes])

        go_file, archA, goaB = self.args
        self.log.info("Loading GO tree from %s..." % go_file)
        self.go_tree = GOTree.from_obo(go_file)
        goa = self.read_goa_file(goaB, codes)

        if not self.options.source_arch and self.options.both:
            self.error("A source architecture file should be specified in 'both' mode")

        if not self.options.source_arch:
            self._transfer_from_same_file(goa, archA)
        elif not self.options.both:
            self._transfer_from_other_file(goa, archA, self.options.source_arch)
        else:
            self._transfer_from_both(goa, archA, self.options.source_arch)


if __name__ == "__main__":
    sys.exit(TransferFunctionFromDomainArch().run())
