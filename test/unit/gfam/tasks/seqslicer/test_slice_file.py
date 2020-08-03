import unittest

from gfam.tasks.seqslicer.slice_file import SliceFile
from test.fixtures.seqslicer_fixture import SeqslicerFixture


class TestSliceFile(unittest.TestCase):
    def setUp(self):
        fixture = SeqslicerFixture()
        self.slice_file = fixture.get_slice_file()
        self.expected_number_of_slices = fixture.get_number_of_slices()
        self.expected_number_of_proteins = fixture.get_number_of_proteins()

    def test_slice_file_should_recover_the_right_numbers(self):
        # arrange
        slice_file_reader = SliceFile()

        # act
        dict_slices = slice_file_reader.get_slice_file_as_dict(self.slice_file)

        # assert
        self.assertEqual(len(dict_slices), self.expected_number_of_proteins)
        self.assertEqual(
            sum([len(x) for x in dict_slices.values()]), self.expected_number_of_slices)
