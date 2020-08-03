from dataclasses import dataclass


@dataclass(frozen=True)
class ArchitectureAssignment:
    protein_id: str
    protein_length: int
    architecture: str
    residues_covered: int
    architecture_detail: str
    residues_covered_by_novel: int
