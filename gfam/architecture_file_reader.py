"""Classes related to reading protein architecture files in GFam"""

from typing import Iterator, List, Tuple

from gfam.utilities.open_anything import open_anything

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

__all__ = ["ArchitectureFileReader"]


class ArchitectureFileReader:
    """Iterates over architectures in a GFam architecture file, returning
    for each architecture a tuple (architecture, [protein_ids]).
    """

    def __init__(self, filename: str, min_coverage_prop: float = 0.0):
        """Constructor

        Parameters
        ----------
        filename : str
            architecture assignment file to read from
        min_coverage_prop : float, optional
            minimum proportion of sequence covered by an architecture to be
            included in the output, by default 0.0
        """
        self.min_cov = min_coverage_prop
        self._f = filename

    def _read_lines_with_minimum_coverage(self) -> Iterator[Tuple[str, str]]:
        for line in open(self._f):
            fields = line.split("\t")
            prot = fields[0]
            arch = fields[3]
            coverage = float(fields[2]) / float(fields[1])
            if coverage >= self.min_cov:
                yield prot, arch

    def __iter__(self):
        prots: List[str] = []
        previous_arch: str = ""
        for prot, arch in self._read_lines_with_minimum_coverage():
            if previous_arch and previous_arch != arch:
                if prots:
                    yield (previous_arch, prots)
                prots = []
            previous_arch = arch
            prots.append(prot)
        if prots and previous_arch != "":
            yield (previous_arch, prots)
