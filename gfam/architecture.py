"""Classes related to reading protein architecture files in GFam"""

import re
from gfam.utils import open_anything

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

__all__ = ["ArchitectureFileReaderPerArch"]


class ArchitectureFileReaderPerArch(object):
    """Iterates over architectures in a GFam architecture file, returning
    for each architecture a tuple (architecture, [protein_ids]). Note
    that the whole architecture file is loaded into memory before
    iterating
    """
    def __init__(self, filename, min_coverage=0.0):
        self.cov = min_coverage
        self._f = open_anything(filename)

    def __iter__(self):
        """Because the file is sorted there is no need to
        sort it by architecture.
        """
        prots = []
        previous_arch = ""
        for line in self._f:
            fields = line.split("\t")
            prot = fields[0]
            coverage = float(fields[2]) / float(fields[1])
            arch = fields[3]
            if previous_arch and previous_arch != arch:
                if prots:
                    yield (previous_arch, prots)
                prots = []
            previous_arch = arch
            if coverage >= self.cov:
                prots.append(prot)
        if prots and previous_arch != "":
            yield (previous_arch, prots)


class ArchitectureDetailsParser(object):
    """Parser of the details regarding a protein architecture from a GFam
    architecture output file
    """
    def __init__(self):
        self.empty = True
        self.arch = ""
        self.coverage = 0.0
        self.protein = ""

    def add_line(self, line):
        blanks = len(line) - len(line.lstrip())
        if blanks == 0:
            # protein name
            self.protein = line.strip()
        elif blanks >= 3:
            # architecture data
            line = line.strip()
            if line.startswith("Coverage:"):
                (_, cov) = line.split()
                self.coverage = float(cov)
            elif re.match(r"^\d+-\s*\d+:", line):
                # architecture line
                domain = ((line.split(':')[1]).split())[0]
                if line.find('InterPro ID:') != -1:
                    line = line.replace(")", "")
                    if line.find('-->') != -1:
                        ipr = line.split('-->')[1].strip()
                    else:
                        ipr = line.split('InterPro ID:')[1].strip()
                    domain = ipr
                if self.arch:
                    self.arch = self.arch + ";" + domain
                else:
                    self.arch = domain
                self.empty = False

    def data_as_tuple(self):
        """ Returns parsed data as a tuple (protein, architecture, coverage)
        """
        return self.protein, self.arch, self.coverage
