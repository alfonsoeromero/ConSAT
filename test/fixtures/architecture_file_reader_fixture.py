import os
from typing import Dict, Set


class ArchitectureFileReaderFixture:
    def __init__(self):
        current_dir = os.path.dirname(__file__)
        self.data_dir = os.path.join(current_dir, os.pardir, "data")

    def get_architecture_assignment_file(self) -> str:
        return os.path.join(self.data_dir, "find_domain_arch100.txt")

    def get_with_at_least_certain_coverage(self, min_cov: float) ->\
            Dict[str, Set[str]]:
        architecture_file = {}
        for line in open(self.get_architecture_assignment_file()):
            if not line.strip():
                continue
            parts = line.split("\t")
            prot, arch = parts[0], parts[3]
            coverage = float(parts[2]) / float(parts[1])
            if coverage < min_cov:
                continue
            if arch not in architecture_file:
                architecture_file[arch] = set()
            architecture_file[arch].add(prot)
        return architecture_file

    def get(self) -> Dict[str, Set[str]]:
        return self.get_with_at_least_certain_coverage(min_cov=0.0)
