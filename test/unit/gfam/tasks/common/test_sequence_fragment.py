import unittest
from gfam.tasks.common.sequence_fragment import SequenceFragment


class TestSequenceFragment(unittest.TestCase):
    def setUp(self):
        self._sut = SequenceFragment("ABCDEF", 1, 23)

    def test_representation_should_be_correct(self):
        # arrange
        expected_str = "ABCDEF:1-23"

        # act
        repr = str(self._sut)

        # assert
        self.assertEqual(expected_str, repr)

    def test_building_from_string_should_works(self):
        # arrange
        other_fragment = SequenceFragment.from_str("ABCDEF:1-23")

        # act
        other_repr = str(other_fragment)
        this_repr = str(self._sut)

        # assert
        self.assertEqual(other_repr, this_repr)

    def test_equality_should_work(self):
        # arrange
        other_fragment = SequenceFragment.from_str("ABCDEF:1-23")

        # act

        # assert
        self.assertEqual(self._sut, other_fragment)

    def test_overlaps_should_work(self):
        # arrange
        other_different_id = SequenceFragment.from_str("AAAA:1-10")
        other_one_residue_overlap = SequenceFragment.from_str("ABCDEF:23-50")
        other_next_but_no_overlap = SequenceFragment.from_str("ABCDEF:24-50")

        # act

        # assert
        self.assertFalse(self._sut.overlaps(other_different_id))
        self.assertTrue(self._sut.overlaps(other_one_residue_overlap))
        self.assertFalse(self._sut.overlaps(other_next_but_no_overlap))

    def test_overlap_proportion_should_work(self):
        # arrange
        other_different_id = SequenceFragment.from_str("AAAA:1-10")
        other_one_residue_overlap = SequenceFragment.from_str("ABCDEF:23-50")
        other_next_but_no_overlap = SequenceFragment.from_str("ABCDEF:24-50")

        # act

        # assert
        self.assertEqual(self._sut.overlap_proportion(other_different_id), 0.0)
        self.assertEqual(self._sut.overlap_proportion(
            other_next_but_no_overlap), 0.0)
        self.assertEqual(self._sut.overlap_proportion(
            other_one_residue_overlap), 1/50)
        self.assertEqual(
            other_one_residue_overlap.overlap_proportion(self._sut), 1/50)
        self.assertEqual(self._sut.overlap_proportion(self._sut), 1)
