import os
from typing import List

from gfam.assignment import EValueFilter


class AssignmentReaderWithFiltersFixture:
    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")
        self.evalue_filter_string = "HMMPfam=0.001;HMMSmart=0.005;0.007"
        self.evalue_filter = {k: float(v) for k, v in
                              [x.split("=") if "=" in x else
                               ("default", x) for x in
                               self.evalue_filter_string.split(";")]}

    def get_ignored_sources(self) -> List[str]:
        return ["PrfScan"]

    def get_evalue_filter(self) -> EValueFilter:
        return EValueFilter.from_string(self.evalue_filter_string)

    def get_input_assignment_file(self) -> str:
        return os.path.join(self.data_dir, "uniprot_trembl.interpro")

    def _is_acceptable(self, line: str) -> bool:
        parts = line.strip().split("\t")
        try:
            evalue = float(parts[8])
        except (ValueError, TypeError):
            return True

        source = parts[3]
        if source in self.evalue_filter:
            return evalue < self.evalue_filter[source]
        else:
            return evalue < self.evalue_filter["default"]

    def get_expected_assignments(self) -> List[str]:
        return [line for line in open(self.get_input_assignment_file())
                if all([y not in line for y in self.get_ignored_sources()])
                and self._is_acceptable(line)]
