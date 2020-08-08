import unittest
from test.fixtures.architecture_file_reader_fixture import \
    ArchitectureFileReaderFixture
from gfam.architecture_file_reader import ArchitectureFileReader


class TestArchitectureFileReader(unittest.TestCase):
    def setUp(self):
        self.fixture = ArchitectureFileReaderFixture()

    def test_no_output_should_be_obtained_with_coverage_over_1(self):
        # arrange
        architecture_file = self.fixture.get_architecture_assignment_file()
        sut = ArchitectureFileReader(architecture_file, 1.01)
        expected_assignments = {}

        # act
        obtained_assignments = {x: y for x, y in sut}

        # assert
        self.assertDictEqual(expected_assignments, obtained_assignments)

    def test_read_with_coverage_should_only_give_a_subset(self):
        # arrange
        min_coverage = 0.3
        architecture_file = self.fixture.get_architecture_assignment_file()
        sut = ArchitectureFileReader(architecture_file, min_coverage)
        expected_assignments = self.fixture.get_with_at_least_certain_coverage(
            min_coverage)

        # act
        obtained_assignments = {x: y for x, y in sut}

        # assert
        print(len(obtained_assignments), len(expected_assignments))
        self.assertSetEqual(set(obtained_assignments.keys()),
                            set(expected_assignments.keys()))
        for k, obtained in obtained_assignments.items():
            self.assertSetEqual(set(obtained), expected_assignments[k])

    def test_read_output_should_be_the_same_as_file_content(self):
        # arrange
        architecture_file = self.fixture.get_architecture_assignment_file()
        sut = ArchitectureFileReader(architecture_file)
        expected_assignments = self.fixture.get()

        # act
        obtained_assignments = {x: y for x, y in sut}

        # assert
        self.assertSetEqual(set(obtained_assignments.keys()),
                            set(expected_assignments.keys()))
        for k, obtained in obtained_assignments.items():
            self.assertSetEqual(set(obtained), expected_assignments[k])
