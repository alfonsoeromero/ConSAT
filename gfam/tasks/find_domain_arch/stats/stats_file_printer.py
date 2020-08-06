from gfam.tasks.find_domain_arch.stats.architecture_stats import\
    ArchitectureStats
from gfam.utils import redirected


class StatsFilePrinter:
    def _get_stats_text(self, stats: ArchitectureStats) -> str:
        text = f"""Domain architectures
               ====================

               Non-empty {stats.total_archs:,}
               Non-empty (when ignoring novel domains): {
               stats.archs_without_novel:,}

               Sequences
               =========

               Total: {stats.total_sequences:,}"""
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
