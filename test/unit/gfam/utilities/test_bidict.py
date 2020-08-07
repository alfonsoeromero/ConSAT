import unittest

from gfam.utilities.bidict import Bidict


class TestBidict(unittest.TestCase):
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
