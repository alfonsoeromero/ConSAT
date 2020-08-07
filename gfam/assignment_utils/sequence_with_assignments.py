import operator
from dataclasses import replace
from typing import List

from gfam.assignment_utils.assignment import Assignment
from gfam.assignment_utils.assignment_overlap_checker import \
    AssignmentOverlapChecker
from gfam.assignment_utils.overlap_type import OverlapType


class SequenceWithAssignments:
    """Class representing a sequence for which some parts are assigned to
    InterPro domains.

    The class has the following fields:

    - ``name``: the name of the sequence

    - ``length``: the number of amino acids in the sequence

    - ``assignments``: a list of `Assignment` instances that describe
      the domain architecture of the sequence

    - ``architecture``: architecture of the sequence, in a string form
      computed from ``assignments`` using the method ``get_string`` of
      ``TreeRepresentation``.

    - ``architecture_pos``: same as the previous one, but containing
      the positions in which each domain spans between parentheses.
    """

    __slots__ = ("name", "length", "assignments", "architecture",
                 "architecture_pos")

    #: The overlap checker used by this instance. This points to
    #: `AssignmentOverlapChecker` by default.
    overlap_checker = AssignmentOverlapChecker

    #: Acceptable overlap types which we allow in assignments.
    #: By default, we consider complete domain insertions and the
    #: absence of overlaps acceptable.
    acceptable_overlaps = set([OverlapType.NO_OVERLAP, OverlapType.INSERTION])

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.assignments = []

    def __len__(self):
        return self.length

    def assign_(self, start, end, domain, source="Novel", *args, **kwds):
        """Assigns a fragment of this sequence to the given domain.
        `start` and `end` are the starting and ending positions, inclusive.
        `domain_name` is the name of the domain, `source` is the assignment
        source (``Novel`` by default).
        """
        assignment = Assignment(id=self.name, start=start, end=end,
                                interpro_id=None, source=source, domain=domain,
                                evalue=None, length=self.length, comment=None)
        return self.assign(assignment, *args, **kwds)

    def assign_list_without_checking_overlap(
            self, assignments: List[Assignment]) -> None:
        """Assigns a list of assigments in batch to the sequence,
            without checking the overlap

        Parameters
        ----------
        assignments : List[Assignment]
            list of assignments to assign to the sequence
        """
        for assignment in assignments:
            self.assign(assignment, overlap_check=False)

    def assign(self, assignment: Assignment, overlap_check=True,
               interpro=None):
        """Assigns a fragment of this sequence using the given assignment.
        If `overlap_check` is ``False``, we will not check for overlaps or
        conflicts with existing assignments.

        Returns ``True`` if the assignment was added, ``False`` if it
        wasn't due to an overlap conflict.
        """
        if ":SF" in assignment.domain:
            # Remove subfamily info from domain
            # pylint: disable-msg=W0212
            new_domain = assignment.domain[0:assignment.domain.index(":SF")]
            assignment = replace(assignment, domain=new_domain)

        if overlap_check:
            overlap_state = self.overlap_checker.check(self, assignment,
                                                       interpro)
            if overlap_state not in self.acceptable_overlaps:
                return False
        self.assignments.append(assignment)
        return True

    def num_covered(self, sources=None):
        """Returns the number of residues covered by the assignments in the
        sequence.

        `sources` specifies the data sources to be included in the coverage
        calculation. If `None`, all the data sources will be considered;
        otherwise it must be a set containing the accepted sources.
        """
        ok = [0] * self.length
        if sources is None:
            for a in self.assignments:
                ok[a.start:(a.end+1)] = [1] * ((a.end+1)-a.start)
        else:
            if isinstance(sources, str):
                sources = [sources]
            for a in self.assignments:
                if a.source in sources:
                    ok[a.start:(a.end+1)] = [1] * ((a.end+1)-a.start)
        return sum(ok)

    def coverage(self, sources=None):
        """Returns the coverage of the sequence, i.e. the fraction of residues
        covered by at least one assignment.

        `sources` specifies the data sources to be included in the coverage
        calculation. If `None`, all the data sources will be considered;
        otherwise it must be a set containing the accepted sources."""
        return self.num_covered(sources) / float(self.length)

    def data_sources(self):
        """Returns the list of data sources that were used in
        this assignment."""
        return sorted(set(a.source for a in self.assignments))

    def domain_architecture(self, sources=None):
        """Returns the domain architecture of the assignment.

        The domain architecture is a list which contains the IDs of the
        assigned regions (domains) in ascending order of their starting
        positions. If `sources` is ``None``, all data sources will be
        considered; otherwise it must be a set or iterable which specifies
        the data sources to be included in the result.
        """
        sorted_assignments = sorted(self.assignments,
                                    key=operator.attrgetter("start"))
        if sources is None:
            return [a.domain for a in sorted_assignments]
        if isinstance(sources, str):
            sources = [sources]
        return [a.domain for a in sorted_assignments if a.source in sources]

    def resolve_interpro_ids(self, interpro):
        """Calls `Assignment.resolve_interpro_ids` on each assignment of this
        sequence"""
        self.assignments = [assignment.resolve_interpro_ids(interpro)
                            for assignment in self.assignments]

    def unassigned_regions(self):
        """Returns a generator that iterates over the unassigned regions
        of the sequence. Each entry yielded by the generator is a tuple
        containing the start and end positions"""
        ok = [True] * (self.length+1)
        for a in self.assignments:
            ok[a.start:(a.end+1)] = [False] * ((a.end+1)-a.start)
        i = 1
        while i <= self.length:
            while i <= self.length and not ok[i]:
                i += 1
            start = i
            if start == self.length:
                break
            while i <= self.length and ok[i]:
                i += 1
            yield start, i - 1
