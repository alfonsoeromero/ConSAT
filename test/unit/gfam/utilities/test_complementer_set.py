from gfam.utilities.complementer_set import ComplementerSet

import unittest


class TestComplementerSet(unittest.TestCase):
    def setUp(self):
        self.sample_element_list = [1, "", 31531.425, ()]

    def test_any_element_should_belong_to_an_empty_complementer_set(self):
        # arrange
        sut = ComplementerSet()

        # act
        membership = [x in sut for x in self.sample_element_list]

        # assert
        self.assertListEqual(membership, [True] * len(membership))

    def test_any_removed_element_should_not_be_in_complementer_set(self):
        # arrange
        sut = ComplementerSet()
        sut -= set([x for x in self.sample_element_list])

        # act
        membership = [x in sut for x in self.sample_element_list]

        # assert
        self.assertListEqual(membership, [False] * len(membership))

    def test_element_removed_several_times_should_raise_exception(self):
        # arrange
        sut = ComplementerSet()

        # act
        sut.remove(2)

        # assert
        with(self.assertRaises(KeyError)):
            sut.remove(2)

    def test_any_unhashable_value_is_always_contained(self):
        # arrange
        sut = ComplementerSet()

        # act
        membership = set() in sut

        # assert
        self.assertTrue(membership)

    def test_difference_should_result_in_intersection_set_or_break_if_not_set(
            self):
        # arrange
        s1 = set([1, 2, 3])

        # act
        result1 = s1 - ComplementerSet()
        result2 = s1 - ComplementerSet([1, 2, 3, 4])

        # assert
        self.assertSetEqual(result1, set())
        self.assertSetEqual(result2, set([1, 2])])
        with self.assertRaises(NotImplementedError):
            2 - ComplementerSet()
