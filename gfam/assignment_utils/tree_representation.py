import operator
from collections import defaultdict
from itertools import combinations

from gfam.assignment_utils.assignment_overlap_checker import \
    AssignmentOverlapChecker
from gfam.assignment_utils.overlap_type import OverlapType


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
