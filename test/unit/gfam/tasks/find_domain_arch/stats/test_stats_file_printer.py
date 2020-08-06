import unittest
from unittest.mock import MagicMock, patch, call

from gfam.tasks.find_domain_arch.stats.architecture_stats import \
    ArchitectureStats
from gfam.tasks.find_domain_arch.stats.stats_file_printer import \
    StatsFilePrinter


class TestStatsFilePrinter(unittest.TestCase):
    def setUp(self):
        self.sample_arch_stats = ArchitectureStats(
            total_archs=10, total_residues=100, total_sequences=12,
            covered_residues=50, covered_residues_nonnovel=33,
            seqs_with_nonempty_domain_arch=7,
            seqs_with_nonempty_domain_arch_ignore_novel=6,
            seqs_with_nonempty_nonnovel_domain_arch=5)

    @patch("builtins.open", new_callable=MagicMock())
    def test_output_should_be_printed_to_stdout_and_file(self, mock_file):
        # arrange
        sut = StatsFilePrinter()
        file_name = "my_file"
        expected_message = sut._get_stats_text(self.sample_arch_stats)

        # act
        sut.print_stats_file(file_name, self.sample_arch_stats)

        # assert
        mock_file.assert_called_with(file_name, "w")
        self.assertIn(call().write(expected_message), mock_file.mock_calls)
        self.assertIn(call().close(), mock_file.mock_calls)
