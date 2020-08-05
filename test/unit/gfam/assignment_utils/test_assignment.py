import os
import unittest

from gfam.assignment_utils.assignment import Assignment
from gfam.interpro import InterPro


class TestAssignment(unittest.TestCase):
    def setUp(self):
        current_dir = os.path.dirname(__file__)
        interpro_file = os.path.join(current_dir, os.pardir,
                                     os.pardir, os.pardir, "data",
                                     "ParentChildTreeFile.txt")
        self.interpro = InterPro.from_file(interpro_file)

    def test_short_repr_should_produce_string(self):
        # arrange
        a = Assignment("my_prot", 120, 5, 50, "HMMPfam",
                       "PFAM0001")

        # act
        short_repr = a.short_repr()
        expected_repr = "PFAM0001(5-50)"

        # assert
        self.assertEqual(short_repr, expected_repr)

    def test_get_assigned_length_should_properly_calculate_length(self):
        # arrange
        a = Assignment("my_prot", 120, 5, 50, "HMMPfam",
                       "PFAM0001")

        # act
        length = a.get_assigned_length()
        expected_length = 46

        # assert
        self.assertEqual(length, expected_length)

    def test_resolve_interpro_ids_should_produce_different_assignment(self):
        # arrange
        a = Assignment("my_prot", 120, 5, 50, "HMMPfam",
                       "PFAM0001", interpro_id="IPR015760")

        # act
        other_assignment = a.resolve_interpro_ids(self.interpro)

        # assert
        self.assertNotEqual(a, other_assignment)
        self.assertEqual(a.id, other_assignment.id)
        self.assertEqual(a.length, other_assignment.length)
        self.assertEqual(a.start, other_assignment.start)
        self.assertEqual(a.end, other_assignment.end)
        self.assertEqual(a.source, other_assignment.source)
        self.assertNotEqual(a.domain, other_assignment.domain)
        self.assertEqual(a.comment, other_assignment.comment)
        self.assertEqual(a.evalue, other_assignment.evalue)
        self.assertEqual(a.interpro_id,
                         other_assignment.interpro_id)
