#!/usr/bin/env python
"""Application that transfers the function of a set of proteins A (from
a GOA file) whose architecture has been computed, from a set of proteins
B, whose architecture is provided.

The GO terms transferred are printed, for each protein id, in a text
readable format
"""

from collections import defaultdict

import operator
import optparse
import sys

from gfam.architecture import ArchitectureFileReader

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"


class TransferFunctionFromDomainArch(CommandLineApp):
    """\
    Usage: %prog [options] architecture_file_A goa_B architecture_file_B output

    Application that transfers the function of a set of proteins A (from
    a GOA file) whose architecture has been computed, from a set of proteins
    B, whose architecture is provided. Results are printed to an `output` file
    """

    short_name = "function_arch"

    def __init__(self, *args, **kwds):
        super(TransferFunctionFromDomainArch, self).__init__(*args, **kwds)
        self.evidence = dict("EXP": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                     "TAS", "IC"],
                             "ALL_BUT_IEA": ["EXP", "IDA", "IPI", "IMP", "IGI",
                                             "IEP", "TAS", "IC", "ISS", "ISO",
                                             "ISA", "ISM", "IGC", "IBA", "IBD",
                                             "IKR", "IRD", "RCA"],
                             "ALL": ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                     "TAS", "IC", "ISS", "ISO", "ISA", "ISM",
                                     "IGC", "IBA", "IBD", "IKR", "IRD", "RCA",
                                     "IEA"]
                             )

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(TransferFunctionFromDomainArch, self).create_parser()
        parser.add_option("-m", "--min-coverage", dest="minimum_coverage",
                          type=float, default=80.0, metavar="PERCENTAGE",
                          config_key="analysis:function_arch/minimum_coverage",
                          help="Minimum % of coverage allowed for an " +
                               "architecture to transfer/receive function")

        parser.add_option("-e", "--ev_codes", dest="ev_codes", type=str,
                          default="EXP", metavar="EXP/ALL_BUT_IEA/ALL",
                          config_key="analysis:function_arch/ev_codes",
                          help="Evidence codes to use: EXP (experimental), " +
                               "ALL_BUT_IEA (all except electronic, IEA," +
                               "NAS and ND), ALL (all evidence codes)")

        return parser

    def read_goa_file(self, goa_file, ev_codes):
        """Reads the GOA file, return a defaultdict that,
        for each protein id it has a set of strings with the
        associated GO terms
        """
        d = defaultdict(set)

        with open(goa_file, "r") as goa:
            for line in goa:
                if not line.startswith("#"):
                    # split line, obtain protein_id, go_term, and ev_code
                    fields = line.split("\t")
                    prot_id, goterm, evcode = fields[1], fields[4], fields[6]
                    if evcode in ev_codes:
                        d[prot_id].add(goterm)
        return d

    def transfer_from_same_file(self, goa, arch_file):
        """ Transfer function from architecture file
        """
        proteins_per_arch = defaultdict(set)
        arch_per_protein = dict()
        for protein, arch, cov in ArchitectureFileReader(arch_file):
            if cov >= self.options.minimum_coverage:
                proteins_per_arch[arch].add(protein)
                arch_per_protein[protein] = arch

        function = defaultdict(list)

        for protein in goa:
            if protein in arch_per_protein:
                arch = arch_per_protein[protein]
                if arch in proteins_per_arch:
                    gos = goa[protein]
                    for target in proteins_per_arch[arch]:
                        if target != protein:
                            for go in gos:
                                function[target].append((go, protein))
        return function

    def transfer_from_other_file(goa, arch_target, arch_source):
        proteins_per_arch = defaultdict(set)

        for protein, arch, cov in ArchitectureFileReader(arch_source):
            if (cov >= self.options.minimum_coverage
               and protein in goa):
                proteins_per_arch[arch].add(protein)

        function = defaultdict(list)

        for protein, arch, cov in ArchitectureFileReader(arch_target):
            if (cov >= self.options.minimum_coverage and
               arch in proteins_per_arch):
                for source in proteins_per_arch[arch]:
                    for go in goa[source]:
                        function[protein].append((go, source))

        return function

    def run_real(self):
        """Runs the applications"""
        if len(self.args) != 4:
            self.error("exactly four input files are expected")

        if self.options.ev_codes not in self.evidence:
            self.error("The three valid types of evidence codes are: " +
                       "EXP, ALL_BUT_IEA, and ALL")

        codes = set(self.evidence[self.options.ev_codes])

        archA, goaB, archB, output = self.args

        goa = self.read_goa_file(goaB, codes)

        if archA == archB:
            function = self.transfer_from_same_file(goa, archA)
        else:
            function = self.transfer_from_other_file(goa, archA, archB)

        with open(output, "w") as out:
            for protein in function:
                out.write("%s\t" % protein)
                for go, source in function[protein]:
                    out.write("%s:%s " % (go, source))

                out.write("\n")


if __name__ == "__main__":
    sys.exit(TransferFunctionFromDomainArch().run())
