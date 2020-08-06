"""Classes related to reading protein architecture files in GFam"""

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
