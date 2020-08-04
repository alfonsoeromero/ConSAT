import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.assignment_source_filter_fixture import \
    AssignmentSourceFilterFixture
from typing import List

from gfam.scripts.assignment_source_filter import AssignmentSourceFilterApp


class TestAssignmentSourceFilter(unittest.TestCase):
    def setUp(self):
        self.fixture = AssignmentSourceFilterFixture()

    def _get_obtained_output(self, sut: AssignmentSourceFilterApp,
                             args: List[str]) -> str:
        """Runs the given `AssignmentSourceFilterApp` `sut` with the
        provided list of arguments `args` returning the output
        as a string"""
        output = io.StringIO()
        # obtain output from find_unassigned and print it to a
        # string stream
        with redirect_stdout(output):
            sut.run(args=args)
        return output.getvalue()

    def test_output_should_match_expected(self):
        # arrange
        args = self.fixture.get_default_args_for_app()
        expected_output = self.fixture.get_expected_assignment_file_as_string()
        _sut = AssignmentSourceFilterApp()

        # act
        obtained_output = self._get_obtained_output(_sut, args)

        # assert
        self.assertEqual(obtained_output, expected_output)
