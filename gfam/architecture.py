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


#class ArchitectureTableFileReaderPerArch(object):
#    def __init__


class ArchitectureFileReaderPerArch(object):
    """Iterates over architectures in a GFam architecture file, returning
    for each architecture a tuple (architecture, [protein_ids]). Note
    that the whole architecture file is loaded into memory before
    iterating
    """
    def __init__(self, filename, min_coverage):
        self.archs = defaultdict(list)
        for protein, arch, cov in ArchitectureFileReader(filename):
            if cov >= min_coverage:
                self.archs[arch].append(protein)

    def __iter__(self):
        for arch in sorted(self.archs):
            yield (arch, self.archs[arch])


class ArchitectureFileReader(object):
    """Iterates over proteins in a GFam architecture file, returning all
    the details (protein, architecture as a string, coverage).
    """

    def __init__(self, filename):
        self._fp = open_anything(filename)

    def architectures(self):
        """A generator that yields tuples (protein, architecture, coverage)
        from a GFam architecture file. Both the `protein` an the `architecture`
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


class ArchitectureParser(object):
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


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "USE: ", sys.argv[0], " architecture_file output_file"
        sys.exit(-1)

    cnt = 0

    with open(sys.argv[2], 'w') as out:
        for protein, arch, cov in ArchitectureFileReader(sys.argv[1]):
            out.write("prot: %s, arch: %s, cov: %f\n" % (protein, arch, cov))
            cnt += 1
        out.write("-------------------------------------------------\n")

    print "File processed, ", cnt, " architectures found."
