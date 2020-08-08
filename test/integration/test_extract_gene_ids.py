import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.extract_gene_ids_fixture import ExtractGeneIdsFixture
from typing import List, Set

from gfam.scripts.extract_gene_ids import ExtractGeneIDsApp


class TestExtractGeneIds(unittest.TestCase):
    def setUp(self):
        self.fixture = ExtractGeneIdsFixture()

    def _get_obtained_output_as_set(self, sut: ExtractGeneIDsApp,
                                    args: List[str]) -> Set[str]:
        """Runs the given `ExtractGeneIDsApp` `sut` with the
        provided list of arguments `args` returning the output
        as a Set of string"""
        output = io.StringIO()
        # obtain output from find_unassigned and print it to a
        # string stream
        with redirect_stdout(output):
            sut.run(args=args)

        # rewind to allow reading it from the beginning
        output_str = output.getvalue()
        return set([x.strip() for x in output_str.split("\n")
                    if x.strip()])

    def test_output_should_match_expected(self):
        # arrange
        args = self.fixture.get_default_args_for_app()
        expected_output_as_set = self.fixture.get_expected_output_as_set()
        _sut = ExtractGeneIDsApp()

        # act
        obtained_output_as_set = self._get_obtained_output_as_set(_sut, args)

        # assert
        self.assertSetEqual(expected_output_as_set,
                            obtained_output_as_set)
