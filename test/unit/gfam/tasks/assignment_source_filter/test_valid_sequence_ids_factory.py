import logging
import unittest
from test.fixtures.assignment_source_filter_fixture import \
    AssignmentSourceFilterFixture

from gfam.tasks.assignment_source_filter.valid_sequence_ids_factory import \
    ValidSequenceIdsFactory
from gfam.utils import complementerset


class TestValidSequenceIdsFactory(unittest.TestCase):
    def setUp(self):
        fixture = AssignmentSourceFilterFixture()
        self.logging = logging.getLogger(__name__)
        self.ids_file = fixture.get_identifiers_file()
        self.expected_ids = {line.strip() for line in open(self.ids_file)}

    def test_get_without_file_should_give_a_complementer_set(self):
        # arrange
        sut = ValidSequenceIdsFactory(self.logging)

        # act
        ids = sut.get()

        # assert
        self.assertIsInstance(ids, complementerset)

    def test_get_with_file_should_give_a_set_of_file_lines(self):
        # arrange
        sut = ValidSequenceIdsFactory(self.logging, self.ids_file)

        # act
        ids = sut.get()

        # assert
        self.assertSetEqual(ids, self.expected_ids)
