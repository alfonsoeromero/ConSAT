import unittest

from gfam.assignment_utils.assignment import Assignment
from gfam.assignment_utils.sequence_with_assignments import \
    SequenceWithAssignments


class TestSequenceWithAssignments(unittest.TestCase):
    def setUp(self):
        pass

    def test_empty_sequence_should_given_one_unassigned_region(self):
        # arrange
        name = "prot"
        length = 100
        sut = SequenceWithAssignments(name, length)

        # act
        regions = list(sut.unassigned_regions())

        # assert
        self.assertEqual(len(regions), 1)
        self.assertEqual(tuple(regions[0]), (1, length))

    def test_len_on_empty_sequence_should_work(self):
        # arrange
        name = "prot"
        length = 100
        sut = SequenceWithAssignments(name, length)

        # act
        found_length = len(sut)

        # assert
        self.assertEqual(length, found_length)

    def test_no_unassigned_regions_should_be_returned_if_sequence_is_covered(
            self):
        # arrange
        name = "prot"
        id = "prot_assignment"
        length = 100
        sut = SequenceWithAssignments(name, length)
        assignment = Assignment(id, 100, 0, 99, "HMMPFam", "domain")

        # act
        sut.assign(assignment)
        unassigned_regions = list(sut.unassigned_regions())

        # assert
        self.assertEqual(len(unassigned_regions), 0)
