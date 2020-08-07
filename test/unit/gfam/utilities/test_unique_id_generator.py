import itertools
import random
import unittest

from gfam.utilities.unique_id_generator import UniqueIdGenerator


class TestUniqueIDGenerator(unittest.TestCase):
    def setUp(self):
        random.seed(42)

    def test_generated_id_should_be_coherent(self):
        # arrange
        sut = UniqueIdGenerator()
        times = 5
        alphabet = "abcdefghijklmopqrstuvwxyz"
        objects = list(alphabet) * times

        # act
        ids = [sut[n] for n in objects]
        mapping = {obj: set(id for _, id in grouped_ids)
                   for obj, grouped_ids in itertools.groupby(zip(objects, ids),
                                                             lambda x: x[0])}

        # assert
        self.assertEqual(len(set(ids)), len(alphabet))
        self.assertSetEqual({len(m) for m in mapping.values()}, set([1]))

    def test_generated_id_should_start_from_given_number(self):
        # arrange
        start_value = 5
        sut = UniqueIdGenerator(id_generator=start_value)
        objects = list("abcde")

        # act
        min_id = min(sut[o] for o in objects)

        # assert
        self.assertEqual(min_id, start_value)

    def test_len_should_give_number_of_unique_generated_ids(self):
        # arrange
        ids_to_generate = 1000
        sut = UniqueIdGenerator()
        objects = [random.randint(1, 50) for _ in range(ids_to_generate)]

        # act
        [sut[o] for o in objects]

        # assert
        self.assertEqual(len(set(objects)), len(sut))

    def test_values_should_give_number_of_mapped_values(self):
        # arrange
        ids_to_generate = 1000
        sut = UniqueIdGenerator()
        objects = [random.randint(1, 50) for _ in range(ids_to_generate)]

        # act
        [sut[o] for o in objects]

        # assert
        self.assertSetEqual(set(sut.values()), set(objects))

    def test_when_generator_is_specified_ids_should_come_from_it(self):
        # arrange
        def _generator() -> int:
            start = 0
            while True:
                yield start
                start += 2
        sut = UniqueIdGenerator(id_generator=_generator())
        ids_to_generate = 500
        objects = [random.randint(1, 50) for _ in range(ids_to_generate)]

        # act
        mapped_objects = [sut[o] for o in objects]

        # assert
        self.assertTrue(all([x % 2 == 0 for x in mapped_objects]))
