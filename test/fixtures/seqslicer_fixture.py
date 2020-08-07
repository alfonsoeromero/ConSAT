import os
from typing import List


class SeqslicerFixture:
    """Fixture to provide better access to testing files, making routes/
        file name independent to the test class"""

    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")

    def get_fasta_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_sprot100.fasta")

    def get_slice_file(self) -> str:
        return os.path.join(self.data_dir, "find_unassigned100.txt")

    def _get_expected_output_file(self) -> str:
        return os.path.join(self.data_dir, "seqslicer.fasta")

    def get_expected_output_as_str(self) -> str:
        with open(self._get_expected_output_file(), "r") as f_in:
            file_content = f_in.read()
        return file_content

    def get_number_of_slices(self) -> int:
        with open(self.get_slice_file()) as f_in:
            num_slices = sum([1 for _ in f_in])
        return num_slices

    def get_number_of_proteins(self) -> int:
        with open(self.get_slice_file()) as f_in:
            num_prots = len(set([f.split()[0]
                                 for f in f_in]))
        return num_prots

    def get_default_args_for_app(self) -> List[str]:
        slice_file = self.get_slice_file()
        fasta_file = self.get_fasta_file()

        args = [slice_file,
                fasta_file,
                "--seq-id-regexp", r"(\w+\|)(?P<id>\w+)(\|\w+)+"]
        return args
