import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.assignment_source_filter_fixture import \
    AssignmentSourceFilterFixture
from typing import List, Union

import pandas as pd
from gfam.scripts.assignment_source_filter import AssignmentSourceFilterApp
from pandas.testing import assert_frame_equal


class TestAssignmentSourceFilter(unittest.TestCase):
    def setUp(self):
        self.fixture = AssignmentSourceFilterFixture()

    def _get_obtained_output(self, sut: AssignmentSourceFilterApp,
                             args: List[str],
                             as_string: bool = True) ->\
            Union[str, pd.DataFrame]:
        """Runs the given `AssignmentSourceFilterApp` `sut` with the
        provided list of arguments `args` returning the output
        as a string or a dataframe, depending on the flag"""
        output = io.StringIO()
        # obtain output from find_unassigned and print it to a
        # string stream
        with redirect_stdout(output):
            sut.run(args=args)

        if as_string:
            return output.getvalue()
        else:
            output.seek(0)
            df_obtained = pd.read_csv(output, sep="\t", header=None)
            return df_obtained

    def test_output_should_match_expected_as_text(self):
        # arrange
        args = self.fixture.get_default_args_for_app()
        expected_output = self.fixture.get_expected_assignment_file_as_string()
        _sut = AssignmentSourceFilterApp()

        # act
        obtained_output = self._get_obtained_output(_sut, args)

        # assert
        self.assertEqual(obtained_output, expected_output)

    def test_output_should_match_expected_as_dataframe(self):
        # arrange
        args = self.fixture.get_default_args_for_app()
        df_expected_output = self.fixture.\
            get_expected_assignment_file_as_dataframe()
        _sut = AssignmentSourceFilterApp()

        # act
        df_obtained_output: pd.DataFrame = self._get_obtained_output(
            _sut, args, as_string=False)

        # assert
        assert_frame_equal(df_obtained_output, df_expected_output,
                           check_names=False, check_like=True)
