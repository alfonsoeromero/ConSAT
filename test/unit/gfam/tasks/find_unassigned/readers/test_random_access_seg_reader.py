import unittest
from collections import defaultdict
from test.fixtures.seg_reader_fixture import SEGReaderFixture
from typing import List

from gfam.tasks.find_unassigned.readers.random_access_seg_reader import \
    RandomAccessSEGReader


class TestRandomAccessSEGReader(unittest.TestCase):
    def setUp(self):
        self.fixture = SEGReaderFixture()
        sequence_id_regexp = r"(\w+\|)(?P<id>\w+)(\|\w+)+"
        self._sut = RandomAccessSEGReader(
            self.fixture.seg_file,
            sequence_id_regexp=sequence_id_regexp)

    def _get_segs_per_protein_using_reader(self, proteins: List[str]) ->\
            defaultdict:
        low_complexity_regions = defaultdict(set)
        for protein in proteins:
            low_complexity_regions[protein] = set(
                self._sut.get_intervals_for_protein(protein))

        return low_complexity_regions

    def test_output_should_match(self):
        # arrange
        expected_output = self.fixture.get_dict_low_complexity_regions()
        proteins = expected_output.keys()

        # act
        obtained_output = self._get_segs_per_protein_using_reader(proteins)

        # assert
        self.assertDictEqual(expected_output, obtained_output)

    def test_non_existing_protein_returns_nothing(self):
        # arrange
        invented_protein = "ABABAB"

        # act
        intervals = list(self._sut.get_intervals_for_protein(
            invented_protein))

        # assert
        self.assertEqual(intervals, [])

    def tearDown(self):
        self._sut.delete_temp_file()
