from typing import Iterator, Set, Tuple

from gfam.assignment_utils.assignment import Assignment
from gfam.assignment_utils.evalue_filter import EValueFilter
from gfam.interpro import AssignmentReader


class AssignmentReaderWithFilters:
    """Similar to `AssignmentReader` but including
        an EvalueFilter and a set of ignored sources.
    """

    def __init__(self, f_name: str, evalue_filter: EValueFilter,
                 ignored_sources: Set[str]):
        """Constructor

        Parameters
        ----------
        f_name : str
            InterPro assignment's file
        evalue_filter : EValueFilter
            Evalue filter
        ignored_sources : Set[str]
            Set with ignored sources
        """
        self.ignored = ignored_sources
        self.evalue_filter = evalue_filter
        self.reader = AssignmentReader(f_name)

    def assignments_and_lines(self) -> Iterator[Tuple[Assignment, str]]:
        """Iterates over an assignment (Interpro) file, producing pairs
        of Assignment objects and lines, that pass the filters defined.

        Returns
        -------
        Tuple[Assignment, str]
            Pair of (Assignment, line) that correspond to an assignment
                passign the above-defined filters

        Yields
        -------
        Iterator[Tuple[Assignment, str]]
            Iterator of the defined pairs
        """
        for assignment, line in self.reader.assignments_and_lines():
            if (assignment.source not in self.ignored) and\
                (assignment.evalue is None or
                 self.evalue_filter.is_acceptable(assignment)):
                yield assignment, line
