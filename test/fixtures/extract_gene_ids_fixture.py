import os
from typing import List, Set


class ExtractGeneIdsFixture:
    """Fixture to provide better access to testing files, making routes/
        file name independent to the test class"""

    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")

    def get_fasta_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_sprot100.fasta")

    def get_ids_file(self) -> str:
        return os.path.join(self.data_dir, "first100IDS.txt")

    def get_expected_output_as_set(self) -> Set[str]:
        return set([x.strip() for x in open(self.get_ids_file()) if x])

    def get_default_args_for_app(self) -> List[str]:
        fasta_file = self.get_fasta_file()

        args = [fasta_file,
                "--seq-id-regexp", r"(\w+\|)(?P<id>\w+)(\|\w+)+"]
        return args
