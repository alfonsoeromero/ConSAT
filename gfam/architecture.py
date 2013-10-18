"""Classes related to reading protein architecture files in GFam"""

import re
import sys
from collections import defaultdict
from gfam.utils import open_anything

try:
    from collections import Mapping
except ImportError:
    from gfam.compat import Mapping

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

__all__ = ["ArchitectureFileReader, ArchitectureFileReaderPerArch"]

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
            prot, coverage, arch = fields[0], float(fields[2])/float(fields[1]), fields[3]
            if previous_arch != "" and previous_arch != arch:
                if prots:
                    yield (previous_arch, prots)
                prots = []
            previous_arch = arch
            if coverage >= self.cov:
                prots.append(prot)
        if prots and previous_arch != "":
            yield (previous_arch, prots)           


class ArchitectureFileReader_old(object):
    """Iterates over proteins in a GFam architecture file, returning all
    the details (protein, architecture as a string, coverage).
    """

    def __init__(self, filename):
        self._fp = open_anything(filename)

    def architectures(self):
        """A generator that yields tuples (protein, architecture, coverage)
        from a GFam architecture file table. Both the `protein` an the `architecture`
        will be strings, and the `coverage` a floating point number between 0.0
        and 100.0
        """
        current_arch = ArchitectureParser()
        for line in self._fp:
            if len(line.strip()) > 0:
                current_arch.add_line(line)
            else:
                if not current_arch.empty:
                    yield current_arch.data_as_tuple()
                current_arch = ArchitectureParser()
        if not current_arch.empty:
            yield current_arch.data_as_tuple()

    def __iter__(self):
        return self.architectures()


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
            elif re.match("^\d+-\s*\d+:", line):
                # architecture line
                domain = ((line.split(':')[1]).split())[0]
                if line.find('InterPro ID:') != -1:
                    line = line.replace(")", "")
                    if line.find('-->') != -1:
                        ipr = line.split('-->')[1].strip()
                    else:
                        ipr = line.split('InterPro ID:')[1].strip()
                    domain = ipr
                if len(self.arch) > 0:
                    self.arch = self.arch + ";" + domain
                else:
                    self.arch = domain
                self.empty = False

    def data_as_tuple(self):
        """ Returns parsed data as a tuple (protein, architecture, coverage)
        """
        return self.protein, self.arch, self.coverage

