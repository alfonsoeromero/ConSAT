import unittest

from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet


class TestSequenceFragmentSet(unittest.TestCase):

    def setUp(self):
        fragments_line = "ABCDEF:1-23  BAADD:40-58   PROTEIN:1-117 PROTEIN:4-24"
        self._sut = SequenceFragmentSet.from_str(fragments_line)

    def test_sizes_should_be_correct_with_repeated_ids(self):
        # arrange
        # act
        num_fragments = self._sut.size()
        num_different_sequences = self._sut.num_different_sequences()

        # assert
        self.assertEqual(num_fragments, 4)
        self.assertEqual(num_different_sequences, 3)

    def test_fragment_set_from_blank_string_should_have_zero_size(self):
        # arrange
        line = ""

        # act
        sut = SequenceFragmentSet(line)

        # assert
        self.assertEqual(sut.size(), 0)
        self.assertEqual(sut.num_different_sequences(), 0)
        self.assertEqual([x for x in sut], [])

    def test_iterator_should_give_same_list(self):
        # arrange
        original_list = list(self._sut.clusters)

        # act
        built_list = [x for x in self._sut]

        # assert
        self.assertListEqual(original_list, built_list)
