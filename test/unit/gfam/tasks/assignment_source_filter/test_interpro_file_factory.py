import unittest
from test.fixtures.assignment_source_filter_fixture import \
    AssignmentSourceFilterFixture

from gfam.tasks.assignment_source_filter.interpro_file_factory import \
    InterproFileFactory


class TestInterproFileFactory(unittest.TestCase):
    def setUp(self):
        self.fixture = AssignmentSourceFilterFixture()

    def test_creating_intepro_without_file_should_work(self):
        # arrange
        # act
        _sut = InterproFileFactory.get_from_file(None)

        # assert
        self.assertEqual(len(_sut.mapping), 0)
        self.assertEqual(len(_sut.tree), 0)

    def test_creating_interpro_from_invalid_file_should_raise_exception(self):
        # arrange
        invalid_file_name = "invalid_file_name"

        # act / assert
        with self.assertRaises(FileNotFoundError):
            InterproFileFactory.get_from_file(invalid_file_name)

    def test_creating_interpro_from_valid_file_name_should_populate_it(self):
        # arrange
        valid_file_name = self.fixture.get_parent_child_interpro_file()

        # act
        _sut = InterproFileFactory.get_from_file(valid_file_name)

        # assert
        self.assertGreater(len(_sut.tree), 0)
        self.assertGreater(len(_sut.mapping), 0)
