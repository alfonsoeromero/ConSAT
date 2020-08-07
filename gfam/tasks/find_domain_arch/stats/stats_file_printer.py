from gfam.tasks.find_domain_arch.stats.architecture_stats import\
    ArchitectureStats
from gfam.utilities.file_utils import redirected


class StatsFilePrinter:
    def _get_stats_text(self, stats: ArchitectureStats) -> str:
        text = f"""Domain architectures
               ====================

               Non-empty {stats.total_archs:,}
               Non-empty (when ignoring novel domains): {
               stats.archs_without_novel:,}

               Sequences
               =========

               Total: {stats.total_sequences:,}
               With at least one domain: {
               stats.seqs_with_nonempty_domain_arch:,} ({
               stats.get_perc_seq_nonempty_arch():.4f}%)
               With at least one non-novel domain: {
               stats.seqs_with_nonempty_domain_arch_ignore_novel:,} ({
               stats.get_perc_seq_nonempty_arch_ignore_novel():.4f}%)
               With at least one domain and no novel domains: {
               stats.seqs_with_nonempty_nonnovel_domain_arch} ({
               stats.get_perc_seq_nonempty_nonnovel_arch():.4f}%)

               Residues
               ========

               Total: {stats.total_residues:,}
               Covered: {stats.covered_residues:,} ({
               stats.get_perc_covered_residues():.4f}%)
               Covered by non-novel: {stats.covered_residues_nonnovel:,} ({
                   stats.get_perc_covered_residues_nonnovel():.4f}%)"""
        return text

    def print_stats_file(self, stats_file_name: str,
                         stats: ArchitectureStats) -> None:
        """Prints the info corresponding to a given stats object

        Parameters
        ----------
        stats_file_name : str
            stats file name where the info will also be written
        stats : ArchitectureStats
            stats object
        """
        stats_file = open(stats_file_name, "w")
        with(redirected(stdout=stats_file)):
            print(self._get_stats_text(stats))
        stats_file.close()
