import os
from typing import List

import pandas as pd


class FindDomainArchFixture:
    """Fixture to provide better access to testing files, making routes/
        file name independent to the test class"""

    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")
        self.columns = ["sequence", "length", "residue s_covered"
                        "architecture", "family_length",
                        "arch_str_pos", "arch_description"]

    def _get_assignment_file(self) -> str:
        return os.path.join(self.data_dir, "assignment_source_filter100.txt")

    def _get_fasta_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_sprot100.fasta")

    def _get_expected_output_file(self) -> str:
        return os.path.join(self.data_dir, "find_domain_arch100.txt")

    def _get_interpro_parent_child_file(self) -> str:
        return os.path.join(self.data_dir, "ParentChildTreeFile.txt")

    def _get_interpro_names_file(self) -> str:
        return os.path.join(self.data_dir, "names.dat.gz")

    def _get_clustering_file(self) -> str:
        return os.path.join(self.data_dir, "cca.txt")

    def get_expected_output_df(self) -> pd.DataFrame:
        df = pd.read_csv(self._get_expected_output_file(),
                         sep="\t",
                         names=self.columns)\
            .sort_values(by=self.columns)\
            .reset_index(drop=True)
        return df

    def get_default_args_for_app(self) -> List[str]:
        assignments_file = self._get_assignment_file()
        fasta_file = self._get_fasta_file()

        interpro_parent_child_file = self._get_interpro_parent_child_file()
        interpro_names_file = self._get_interpro_names_file()
        clustering_file = self._get_clustering_file()

        args = [assignments_file,
                clustering_file,
                "-s", "4",
                "-S", fasta_file,
                "-i", interpro_parent_child_file,
                "-n", interpro_names_file,
                "--max-overlap", "20",
                "--seq-id-regexp", r"(\w+\|)(?P<id>\w+)(\|\w+)+"]
        return args
