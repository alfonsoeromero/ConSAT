import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.seqslicer_fixture import SeqslicerFixture
from typing import List

from gfam.scripts.seqslicer import SeqSlicerApp


class TestSeqSlicer(unittest.TestCase):
    def setUp(self):
        self.fixture = SeqslicerFixture()

    def _get_obtained_output_as_str(self, sut: SeqSlicerApp,
                                    args: List[str]) -> str:
        """Runs the given `SeqSlicerApp` `sut` with the
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
        expected_str = self.fixture.get_expected_output_as_str()
        _sut = SeqSlicerApp()

        # act
        obtained_str = self._get_obtained_output_as_str(_sut, args)

        # assert
        self.assertEqual(obtained_str, expected_str)
