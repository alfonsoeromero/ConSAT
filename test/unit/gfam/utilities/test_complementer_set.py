import unittest

from gfam.utilities.complementer_set import ComplementerSet


class TestComplementerSet(unittest.TestCase):
    def test_contains_should_work(self):
        # arrange
        sut = ComplementerSet([1, 2])

        # act
        v0: bool = (sut in sut)
        v1: bool = ("abc" in sut)
        v2: bool = (1 in sut)

        # assert
        self.assertTrue(v0)
        self.assertTrue(v1)
        self.assertFalse(v2)

    def test_equality_operator_should_work(self):
        # arrange
        sut = ComplementerSet([1, 2])

        # act
        v0: bool = (sut == ComplementerSet([1, 2]))
        v1: bool = (sut == ComplementerSet([1, 2, 3]))
        v2: bool = (sut == 1)

        # assert
        self.assertTrue(v0)
        self.assertFalse(v1)
        self.assertFalse(v2)

    def test_difference_update_should_work(self):
        # arrange
        sut = ComplementerSet([1, 2])

        # act
        sut.difference_update([4, 5])
        sut.difference_update([2], [1, 6], [7, 5, "spam"])

        # assert
        self.assertFalse(any(item in sut
                             for item in [2, 1, 6, 7, 5, "spam", 4]))

    def test_discard_should_work(self):
        # arrange
        sut = ComplementerSet()

        # act
        sut.discard(2)

        # assert
        self.assertEqual(sut, ComplementerSet([2]))
