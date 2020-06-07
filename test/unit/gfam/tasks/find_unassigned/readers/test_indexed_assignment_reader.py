import unittest
from test.fixtures.find_unassigned_fixtures import FindUnassignedFixture
from typing import Dict, List, Set

from gfam.tasks.find_unassigned.readers.\
    indexed_assignment_reader import IndexedAssignmentReader


class TestIndexedAssignmentReader(unittest.TestCase):
    def setUp(self):
        self.fixture = FindUnassignedFixture()
        self.assignment_file = self.fixture.get_assignment_file()
        self.sut = IndexedAssignmentReader(self.assignment_file)

    def _get_lines_per_protein_expected(self) -> Dict[str, Set[str]]:
        d = {}
        for line in open(self.assignment_file):
            protein_id = line.split()[0]
            if protein_id not in d:
                d[protein_id] = set()
            d[protein_id].add(line)
        return d

    def _get_lines_per_protein_reader(self,
                                      reader: IndexedAssignmentReader,
                                      protein_ids: List[str]) ->\
            Dict[str, Set[str]]:
        d = {}
        for protein_id in protein_ids:
            d[protein_id] = set(
                reader.get_lines_with_assignments_for_protein(protein_id))
        return d

    def test_non_existing_protein_returns_nothing(self):
        # arrange
        invented_protein = "ABABAB"

        # act
        lines = list(
            self.sut.get_lines_with_assignments_for_protein(invented_protein))

        # assert
        self.assertEqual(lines, [])

    def test_output_should_match_expected(self):
        # arrange
        expected_output = self._get_lines_per_protein_expected()

        # act
        obtained_output = self._get_lines_per_protein_reader(
            self.sut, expected_output.keys())

        # assert
        self.assertEqual(len(expected_output), len(obtained_output))
        self.assertDictEqual(expected_output, obtained_output)

    def tearDown(self):
        self.sut.delete_temp_file()
