import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.find_unassigned_fixtures import FindUnassignedFixture
from typing import List

import pandas as pd
from pandas.testing import assert_frame_equal

from gfam.scripts.find_unassigned import FindUnassignedApp


class TestFindUnassigned(unittest.TestCase):
    def setUp(self):
        self.fixture = FindUnassignedFixture()

    def _get_obtained_df(self, sut: FindUnassignedApp,
                         args: List[str]) -> pd.DataFrame:
        """Runs the given `FindUnassignedApp` `sut` with the
        provided list of arguments `args` returning the output
        as a dataframe (sequence, left, right)"""
        output = io.StringIO()
        # obtain output from find_unassigned and print it to a
        # string stream
        with redirect_stdout(output):
            sut.run(args=args)
        # rewind to allow reading it from the beginning
        output.seek(0)
        df_obtained = pd.read_csv(output, sep="\t",
                                  names=self.fixture.columns)\
            .sort_values(by=self.fixture.columns)\
            .reset_index(drop=True)
        return df_obtained

    def test_output_should_match_expected(self):
        # arrange
        args = self.fixture.get_default_args_for_app()
        expected_df = self.fixture.get_expected_output_df()
        _sut = FindUnassignedApp()

        # act
        obtained_df = self._get_obtained_df(_sut, args)

        # assert
        assert_frame_equal(obtained_df, expected_df)
