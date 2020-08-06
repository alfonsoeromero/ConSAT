"""Classes corresponding to domain assignments (`Assignment`) and
sequences with corresponding domain assignments (`SequenceWithAssignments`).
"""
import operator
from collections import defaultdict
from dataclasses import replace
from enum import Enum
from itertools import combinations

from gfam.assignment_utils.assignment import Assignment

__author__ = "Tamas Nepusz"
__email__ = "tamas@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2010, Tamas Nepusz"
__license__ = "GPL"

__all__ = ["AssignmentOverlapChecker", "OverlapType",
           "SequenceWithAssignments", "EValueFilter", "TreeRepresentation"]

# pylint: disable-msg=C0103,E1101
# C0103: invalid name
# E1101: instance has no 'foo' member. Pylint does not know namedtuple
# trickery.


class OverlapType(Enum):
    """Enum describing the different overlap types that can be detected
    by `AssignmentOverlapChecker`. See the documentation of
    `AssignmentOverlapChecker.check_single` for more details.
    """
    NO_OVERLAP: str = "NO_OVERLAP"
    INSERTION: str = "INSERTION"
    INSERTION_DIFFERENT: str = "INSERTION_DIFFERENT"
    SYNONYM_INSERTION: str = "SYNONYM_INSERTION"
    DIFFERENT: str = "DIFFERENT"
    OVERLAP: str = "OVERLAP"


class TreeRepresentation(object):
    """A class for printing a nice tree-representation of the domains
    composing the architecture of a protein. "tree" is a list where
    each member is a tuple (a, b). The first component, a, is an assignment,
    and the second one, b, is a list (possibly empty) containing the
    "descendants" of that tree node
    """

    def __init__(self, assignments, interpro=None):
        self.interpro = interpro
        self.assignments = assignments
        self.__find_parents(self.assignments, interpro)
        self.tree = self.__get_tree_representation()

    def __find_parents(self, assignments, interpro):
        """Finds the parents of the assignment. The `parents` dictionary
        will store all the possible parents of a domain, if any. Note that
        a domain within a nested insertion will have two parents, which would
        be effectively a branch A -> B -> C.
        """
        self.parents = defaultdict(list)
        for ass1, ass2 in combinations(assignments, 2):
            overlap = AssignmentOverlapChecker.check_single(ass1, ass2,
                                                            interpro)
            if overlap == OverlapType.INSERTION:
                if ass1.get_assigned_length() >= ass2.get_assigned_length():
                    self.parents[ass2].append(ass1)
                else:
                    self.parents[ass1].append(ass2)

    def __get_tree_representation(self):
        li = []
        s = sorted(self.assignments, key=operator.attrgetter("start"))
        for ass in sorted(s, key=operator.attrgetter("length"), reverse=True):
            if ass not in self.parents:
                li.append((ass, []))
            else:
                # a recursively insertion, one level per parent
                current = li
                explored = set()
                level = 0
                pars = self.parents[ass]

                while level < len(pars):
                    for parent in pars:
                        cur = [domain for domain, _ in current]
                        if parent not in explored and parent in cur:
                            _, current = current[cur.index(parent)]
                            explored.add(parent)
                            break
                    level += 1
                current.append((ass, []))
        return li

    def get_string_positions(self, li=None):
        """Gets the tree representation as one string listing the positions
        and the models
        """
        if not self.assignments:
            return "NO_ASSIGNMENT"
        if li is None:
            li = self.tree
        if len(li) == 1:
            assignment, children = li[0]
            rep = assignment.short_repr()
            if children:
                rep += "{" + self.get_string_positions(children) + "}"
            return rep
        else:
            return ";".join([self.get_string_positions([ass]) for ass in li])

    def get_string(self, li=None):
        """Gets the tree representation as one string where only the models
        are listed
        """
        if not self.assignments:
            return "NO_ASSIGNMENT"
        if li is None:
            li = self.tree
        if len(li) == 1:
            assignment, children = li[0]
            rep = assignment.domain
            if children:
                rep += "{" + self.get_string(children) + "}"
            return rep
        else:
            return ";".join([self.get_string([ass]) for ass in li])


class AssignmentOverlapChecker(object):
    """Static class that contains the central logic of determining
    whether an assignment can be added to a partially assigned
    `SequenceWithAssignments`.

    The class has a class variable named `max_overlap` which stores
    the maximum allowed overlap size. This is 20 by default.
    """
    #: The maximum allowed overlap size.
    max_overlap = 20

    #: This is the minimum number of residues that a domain with the same
    #: InterPro id can leave uncovered in the parent until it starts.
    #: Otherwise it is discarded and considered as the same domain. Therefore,
    #: if two domains with the same InterPro are inserted and the beginning of
    #: both has a difference between 0 and 20 residues, the insertion is
    #: labeled as a synonym insertion
    min_parent_inserted_size = 20

    #: Here we give a scale of priorities for the different overlap types.
    #: When checking an assignment against all previously added assignments
    #: we will first check if there is overlap of the first kind
    #: (`OverlapType.OVERLAP`) with any of the already added members.
    #: If so, this overlap type is returned. If not, we pass on to the next
    #: `OverlapType`, until all the types have been tested. This is done to
    #: avoid unwanted behaviour. For instance, an assignment C might be found
    #: as an insertion in a large assignment A, which also contains an
    #: assignment B already inserted, overlapping with C. If we checked only
    #: if there were an insertion, A would be accepted, as no further checks
    #: are done.
    priority = [OverlapType.OVERLAP,
                OverlapType.DIFFERENT, OverlapType.SYNONYM_INSERTION,
                OverlapType.INSERTION_DIFFERENT,
                OverlapType.INSERTION, OverlapType.NO_OVERLAP]

    @classmethod
    def check(cls, sequence, assignment, interpro):
        """Checks whether an `assignment` can be added to a partially
        assigned `sequence`. `sequence` must be an instance of
        `SequenceWithAssignments`, `assignment` must be an instance
        of `Assignment`. InterPro is an interpro data structure
        (we use it to check for similar IPR domains).

        The most strict type of overlap (following the order defined
        in "priority" is returned)
        """
        overlaps_found = set(cls.check_single(assignment, other_assignment,
                                              interpro)
                             for other_assignment in sequence.assignments)
        for overlap_type in cls.priority:
            if overlap_type in overlaps_found:
                return overlap_type
        return OverlapType.NO_OVERLAP

    @classmethod
    def similar(cls, ipr1, ipr2, interpro):
        if ipr1 == ipr2:
            return True
        else:
            return interpro.tree.get_most_remote_ancestor(ipr1) ==\
                interpro.tree.get_most_remote_ancestor(ipr2)

    @classmethod
    def check_single(cls, assignment, other_assignment, interpro=None):
        """Checks whether the given `assignment` overlaps with another
        assignment `other_assignment`. Returns one of the following:

        - `OverlapType.NO_OVERLAP`: there is no overlap between the
          two given assignments

        - `OverlapType.SYNONYM_INSERTION`: `assignment` is inserted into
          `other_assignment` or vice versa, but they share the same InterPro
           term, and very few amino acids in the parent (less than
           `min_parent_inserted_size`) are covered in the starting fragment.
           Therefore we might think it is not a real insertion, but two
           expressions of the same domain (belonging for example to two
           different signature databases).

        - `OverlapType.INSERTION`: there is a complete domain insertion in
          either direction

        - `OverlapType.INSERTION_DIFFERENT`: `assignment` is inserted into
          `other_assignment` or vice versa, but they have different data
          sources.

        - `OverlapType.DIFFERENT`: `other_assignment` overlaps with
        `assignment` partially, but they have different data sources

        - `OverlapType.OVERLAP`: `other_assignment` overlaps with `assignment`
          partially, they have the same data source, but the size of
          the overlap is larger than the maximum allowed overlap specified
          in `AssignmentOverlapChecker.max_overlap`.
        """
        start, end = assignment.start, assignment.end
        other_start, other_end = other_assignment.start, other_assignment.end

        if ((other_start <= start and other_end >= end)
                or (other_start >= start and other_end <= end)):
            # Domains are one inside the other: could be INSERTION,
            # INSERTION_DIFFERENT,
            # SYNONYM_INSERTION or OVERLAP
            similar_iprs = (cls.similar(assignment.interpro_id,
                                        other_assignment.interpro_id,
                                        interpro) or
                            assignment.domain == other_assignment.domain)

            if abs(other_start-start) + abs(other_end-end)\
               < cls.min_parent_inserted_size:
                if similar_iprs and\
                   (abs(other_start - start) < cls.min_parent_inserted_size):
                    return OverlapType.SYNONYM_INSERTION
                else:
                    return OverlapType.OVERLAP
            elif other_assignment.source == assignment.source:
                # This is a valid domain insertion, other_assignment is
                # inserted into assignment
                if not similar_iprs or\
                   (abs(other_start - start) >= cls.min_parent_inserted_size):
                    return OverlapType.INSERTION
                else:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.INSERTION_DIFFERENT

        if (other_start <= start and start <= other_end <= end) or\
           (start <= other_start and other_start <= end <= other_end):
            # Domains could overlap (either A B A' B' or B A B' A')
            if other_assignment.source == assignment.source:
                # This is a partial overlap
                overlap_size = min(other_end - start, end - other_start) + 1
                if overlap_size > cls.max_overlap:
                    return OverlapType.OVERLAP
            else:
                return OverlapType.DIFFERENT

        return OverlapType.NO_OVERLAP


class SequenceWithAssignments(object):
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


class EValueFilter(object):
    """Given an `Assignment`, this filter tells whether the assignment's
    E-value is satisfactory to accept it.

    The filter supports different E-values for different data sources. By
    default, the E-value threshold is infinity for all data sources.
    """

    def __init__(self):
        self.default_e_value = float('inf')
        self.thresholds = {}

    def set_threshold(self, source, evalue):
        """Sets the E-value threshold for the given data source."""
        self.thresholds[source] = evalue

    def is_acceptable(self, assignment):
        """Checks whether the given assignment is acceptable.

        This method looks up the E-value threshold corresponding to the
        ``source`` of the assignment and returns ``True`` if the E-value
        of the assignment is less than the threshold, ``False``
        otherwise."""
        threshold = self.thresholds.get(assignment.source,
                                        self.default_e_value)
        return assignment.evalue <= threshold

    @classmethod
    def from_string(cls, description):
        """Constructs an E-value filter from a string description that
        can be used in command line arguments and configuration files.

        The string description is a semicolon-separated list of
        source-threshold pairs. For instance, the following is a valid
        description giving an E-value of 0.001 for HMMPfam sources,
        0.005 for HMMSmart sources and 0.007 for everything else::

            HMMPfam=0.001;HMMSmart=0.005;0.007

        The last entry denotes the default E-value; in particular, if a
        number is not preceded by a source name, it is assumed to be a
        default E-value. If no default E-value is given, infinity will
        be used.
        """
        result = cls()
        for part in description.split(";"):
            part = part.strip()
            if "=" in part:
                source, evalue = part.split("=", 1)
                evalue = float(evalue)
                result.set_threshold(source, evalue)
            else:
                result.default_e_value = float(part)
        return result
