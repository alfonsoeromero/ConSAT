import unittest
from test.fixtures.assignment_reader_with_filters_fixture import \
    AssignmentReaderWithFiltersFixture
from gfam.tasks.assignment_source_filter.assignment_reader_with_filters import\
    AssignmentReaderWithFilters


class TestAssignmentReaderWithFilters(unittest.TestCase):
    def setUp(self):
        self.fixture = AssignmentReaderWithFiltersFixture()

    def test_assignment_reader_should_filter_with_given_filters(self):
        # arrange
        assignment_file = self.fixture.get_input_assignment_file()
        ignored_sources = self.fixture.get_ignored_sources()
        evalue_filter = self.fixture.get_evalue_filter()
        reader = AssignmentReaderWithFilters(assignment_file,
                                             evalue_filter,
                                             ignored_sources)
        expected_lines = self.fixture.get_expected_assignments()

        # act
        obtained_lines = [line for _, line in reader.assignments_and_lines()]

        # assert
        self.assertEqual(len(expected_lines), len(obtained_lines))
        self.assertListEqual(expected_lines, obtained_lines)
