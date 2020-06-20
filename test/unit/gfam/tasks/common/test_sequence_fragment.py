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
