from dataclasses import dataclass


@dataclass
class ArchitectureStats:
    """Store all the counts made at the end
        of the find_domain_arch apps
    """
    total_archs: int = 0
    total_sequences: int = 0
    total_residues: int = 0

    archs_without_novel: int = 0

    seqs_with_nonempty_domain_arch: int = 0
    seqs_with_nonempty_domain_arch_ignore_novel: int = 0
    seqs_with_nonempty_nonnovel_domain_arch: int = 0

    covered_residues: int = 0
    covered_residues_nonnovel: int = 0

    def get_perc_seq_nonempty_arch(self) -> float:
        return self.seqs_with_nonempty_domain_arch * 100 /\
            self.total_sequences

    def get_perc_seq_nonempty_arch_ignore_novel(self) -> float:
        return self.seqs_with_nonempty_domain_arch_ignore_novel * 100 /\
            self.total_sequences

    def get_perc_seq_nonempty_nonnovel_arch(self) -> float:
        return self.seqs_with_nonempty_nonnovel_domain_arch * 100 /\
            self.total_sequences

    def get_perc_covered_residues(self) -> float:
        return 100 * self.covered_residues / self.total_residues

    def get_perc_covered_residues_nonnovel(self) -> float:
        return 100 * self.covered_residues_nonnovel / self.total_residues
