import os
from typing import List

import pandas as pd


class FindUnassignedFixture:
    """Fixture to provide better access to testing files, making routes/
        file name independent to the test class"""

    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")
        self.columns = ["sequence", "start", "end"]

    def get_assignment_file(self) -> str:
        return os.path.join(self.data_dir, "assignment_source_filter100.txt")

    def _get_low_complexity_regions_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_sprot100.mask")

    def _get_fasta_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_sprot100.fasta")

    def _get_expected_output_file(self) -> str:
        return os.path.join(self.data_dir, "find_unassigned100.txt")

    def get_expected_output_df(self) -> pd.DataFrame:
        df = pd.read_csv(self._get_expected_output_file(),
                         sep="\t",
                         names=self.columns)\
            .sort_values(by=self.columns)\
            .reset_index(drop=True)
        return df

    def get_default_args_for_app(self) -> List[str]:
        assignments_file = self.get_assignment_file()
        fasta_file = self._get_fasta_file()
        seg_file = self._get_low_complexity_regions_file()

        args = [assignments_file,
                "-S", fasta_file,
                "--low-complexity-regions-file", seg_file,
                "--seq-id-regexp", r"(\w+\|)(?P<id>\w+)(\|\w+)+",
                "--max-overlap", "20",
                "--min-fragment-length", "75",
                "--min-length", "30"]
        return args
