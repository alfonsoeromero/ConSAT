import unittest
from gfam.utils import goid_to_num, num_to_goid


class TestGOIdToNum(unittest.TestCase):
    def test_should_return_correct_number(self):
        # arrange
        expected_pairs = {"GO:0001308": 1308,
                          "GO:0000001": 1,
                          "GO:1234567": 1234567}

        # act
        obtained_pairs = {x: goid_to_num(x) for x in expected_pairs}

        # assert
        self.assertDictEqual(expected_pairs, obtained_pairs)


class TestNumToGOID(unittest.TestCase):
    def test_should_return_correct_number(self):
        # arrange
        expected_pairs = {1308: "GO:0001308",
                          1: "GO:0000001",
                          1234567: "GO:1234567"}

        # act
        obtained_pairs = {x: num_to_goid(x) for x in expected_pairs}

        # assert
        self.assertDictEqual(expected_pairs, obtained_pairs)
