import unittest

from gfam.tasks.find_domain_arch.stats.architecture_stats import \
    ArchitectureStats


class TestArchitectureStats(unittest.TestCase):
    def test_architecture_stats_with_totals_1_should_give_0_everywhere_else(
            self):
        # arrange
        sut = ArchitectureStats(
            total_archs=1, total_residues=1, total_sequences=1
        )
        expected_percentages = [0, 0, 0, 0, 0]

        # act
        percentages = [sut.get_perc_covered_residues(),
                       sut.get_perc_covered_residues_nonnovel(),
                       sut.get_perc_seq_nonempty_arch(),
                       sut.get_perc_seq_nonempty_arch_ignore_novel(),
                       sut.get_perc_seq_nonempty_nonnovel_arch()]

        # assert
        self.assertListEqual(percentages, expected_percentages)

    def test_architecture_stats_should_give_correct_percentages(self):
        # arrange
        sut = ArchitectureStats(
            total_archs=10, total_residues=100, total_sequences=12,
            covered_residues=50, covered_residues_nonnovel=33,
            seqs_with_nonempty_domain_arch=7,
            seqs_with_nonempty_domain_arch_ignore_novel=6,
            seqs_with_nonempty_nonnovel_domain_arch=5)
        expected_percentages = [50, 33, 100*7/12, 50, 125/3]

        # act
        percentages = [sut.get_perc_covered_residues(),
                       sut.get_perc_covered_residues_nonnovel(),
                       sut.get_perc_seq_nonempty_arch(),
                       sut.get_perc_seq_nonempty_arch_ignore_novel(),
                       sut.get_perc_seq_nonempty_nonnovel_arch()]

        # assert
        self.assertListEqual(percentages, expected_percentages)
