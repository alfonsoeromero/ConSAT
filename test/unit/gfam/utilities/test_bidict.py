import random
import unittest
from collections import defaultdict

from gfam.utilities.bidict import Bidict


class TestBidict(unittest.TestCase):
    def setUp(self):
        random.seed(42)

    def test_add_left_twice_with_same_left_elem_should_return_set_right(self):
        # arrange
        same_key = "foo"
        different_vals = ["bar", "bzr"]
        sut = Bidict()

        # act
        for val in different_vals:
            sut.add_left(same_key, val)
        elems_right = sut.get_left(same_key)
        elems_left = {val: sut.get_right(val)
                      for val in different_vals}

        # assert
        self.assertSetEqual(elems_right, set(different_vals))
        self.assertDictEqual(elems_left, {val: set([same_key])
                                          for val in different_vals})
        self.assertEqual(sut.len_left(), 1)

    def test_constructing_with_wrong_type_should_raise_exception(self):
        # arrange
        non_iterable = 1234
        # act /  assert
        with(self.assertRaises(TypeError)):
            Bidict(non_iterable)

    def test_constructing_with_dict_should_work(self):
        # arrange
        dictionary = {"a": [0, 1, 2], "b": [3, 4, 5]}
        expected_left = {x: set(y) for x, y in dictionary.items()}
        expected_right = defaultdict(set)
        for x, y in dictionary.items():
            for y_i in y:
                expected_right[y_i].add(x)
        # act
        sut = Bidict(dictionary)

        # assert
        self.assertDictEqual(sut.left, expected_left)
        self.assertDictEqual(sut.right, expected_right)

    def test_iteritems_should_give_all_left_elems(self):
        # arrange
        num_generated = 1000
        left_keys = [random.randint(1, 50) for _ in range(num_generated)]
        right_values = [random.choice(range(5)) for _ in range(num_generated)]

        # act
        sut = Bidict()
        for left, right in zip(left_keys, right_values):
            sut.add_left(left, right)

        # assert
        self.assertListEqual(sorted([x for (x, _) in sut.iteritems_left()]),
                             sorted(set(left_keys)))

    def test_copy_constructed_bidict_should_work_the_same(self):
        same_key = "foo"
        different_vals = ["bar", "bzr"]
        source = Bidict()

        # act
        for val in different_vals:
            source.add_left(same_key, val)
        sut = Bidict(items=source)

        # assert
        self.assertDictEqual(sut.right, source.right)
        self.assertDictEqual(sut.left, source.left)
