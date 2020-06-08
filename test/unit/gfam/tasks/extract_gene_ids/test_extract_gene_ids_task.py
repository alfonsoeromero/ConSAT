import io
import unittest
from contextlib import redirect_stdout
from test.fixtures.find_unassigned_fixtures import FindUnassignedFixture
from typing import List

from gfam.tasks.extract_gene_ids.task import ExtractGeneIDsTask


class TestExtractGeneIDsTask(unittest.TestCase):
    def setUp(self):
        fixture = FindUnassignedFixture()
        self.fasta_file = fixture.get_fasta_file()
        self._sut = ExtractGeneIDsTask()

    def _get_obtained_output(self) -> List[str]:
        output = io.StringIO()
        # obtain output from find_unassigned and print it to a
        # string stream
        with redirect_stdout(output):
            self._sut.process_file(self.fasta_file,
                                   r"(\w+\|)(?P<id>\w+)(\|\w+)+")
        # rewind to allow reading it from the beginning
        output.seek(0)
        obtained_prot_ids = [x.strip() for x in output]
        return obtained_prot_ids

    def test_output_should_be_correct(self):
        # arrange
        expected_output = [x.strip().split("|")[1]
                           for x in open(self.fasta_file)
                           if x and x.startswith(">")]

        # act
        obtained_output = self._get_obtained_output()

        # assert
        self.assertEqual(expected_output, obtained_output)
