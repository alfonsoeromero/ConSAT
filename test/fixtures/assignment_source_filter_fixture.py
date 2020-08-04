import os
from typing import List


class AssignmentSourceFilterFixture:
    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")

    def get_default_args_for_app(self) -> List[str]:
        assignment_file: str = self.get_input_assignment_file()
        interpro_parent_child_file: str = self.get_parent_child_interpro_file()

        args = [assignment_file,
                "-x", "HAMAP PatternScan FPrintScan Seg Coil",
                "-e", "1e-3;superfamily=inf;HMMPanther=inf;Gene3D=inf;HMMPIR=inf",
                "-i", interpro_parent_child_file,
                "--max-overlap", "20"]
        return args

    def get_expected_assignment_file_as_string(self) -> str:
        output_file = os.path.join(
            self.data_dir, "assignment_source_filter100.txt")
        with open(output_file, "r") as f_in:
            out_str = f_in.read()
        return out_str

    def get_input_assignment_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_trembl.interpro")

    def get_identifiers_file(self) -> str:
        return os.path.join(self.data_dir, "first100IDS.txt")

    def get_parent_child_interpro_file(self) -> str:
        return os.path.join(self.data_dir, "ParentChildTreeFile.txt")
