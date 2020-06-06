import unittest
from test.fixtures.find_unassigned_fixtures import FindUnassignedFixture
from typing import Dict, List

from gfam.interpro import Assignment, AssignmentReader
from gfam.tasks.find_unassigned.indexed_readers\
         .random_access_assignment_reader import \
    RandomAccessAssignmentReader


class TestRandomAccessAssignmentReader(unittest.TestCase):
    def setUp(self):
        self.fixture = FindUnassignedFixture()
        self.assignment_file = self.fixture.get_assignment_file()
        self.sut = RandomAccessAssignmentReader(self.assignment_file)

    def _get_assignment_per_protein_manually(self) -> Dict[str, Assignment]:
        ret = {}
        for assignment in AssignmentReader(self.assignment_file):
            prot_id = assignment.id
            if prot_id not in ret:
                ret[prot_id] = set()
            ret[prot_id].add(assignment)
        return ret

    def _get_assignment_per_protein_reader(
            self, reader: RandomAccessAssignmentReader,
            proteins: List[str]) ->\
            Dict[str, Assignment]:
        ret = {}
        for prot in proteins:
            ret[prot] = set(self.sut.assignments_for_protein(prot))
        return ret

    def test_output_should_match(self):
        # arrange
        expected_output = self._get_assignment_per_protein_manually()
        proteins = expected_output.keys()

        # act
        obtained_output = self._get_assignment_per_protein_reader(self.sut,
                                                                  proteins)

        # assert
        self.assertDictEqual(expected_output, obtained_output)

    def test_non_existing_protein_returns_nothing(self):
        # arrange
        invented_protein = "ABABAB"

        # act
        assignments = list(self.sut.assignments_for_protein(invented_protein))

        # assert
        self.assertEqual(assignments, [])

    def tearDown(self):
        self.sut.delete_temp_file()
