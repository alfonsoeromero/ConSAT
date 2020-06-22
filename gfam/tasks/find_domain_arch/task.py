import logging
from typing import Dict, List

from gfam import fasta
from gfam.assignment import AssignmentOverlapChecker
from gfam.interpro import InterPro, InterProNames
from gfam.sequence import SeqRecord
from gfam.tasks.base import LoggedTask
from gfam.tasks.find_unassigned.readers.random_access_assignment_reader import \
    RandomAccessAssignmentReader
from gfam.utils import open_anything, redirected
from gfam.tasks.find_domain_arch.clustering_file import ClusteringFile


class FindDomainArchTask(LoggedTask):
    def __init__(self,
                 max_overlap: int,
                 interpro_names_file: str,
                 min_cluster_size: int,
                 sequence_id_regexp: str = "",
                 interpro_parent_child_file: str = "",
                 logger: logging.Logger = None):
        super().__init__(logger)

        AssignmentOverlapChecker.max_overlap = max_overlap

        if not interpro_parent_child_file:
            self.log.info("Loading InterPro parent-child"
                          " assignments from {interpro_parent_child_file}...")
            self.interpro = InterPro.FromFile(interpro_parent_child_file)
        else:
            self.interpro = InterPro()
        self.interpro_names = InterProNames.FromFile(interpro_names_file)
        self.sequence_id_regexp = sequence_id_regexp
        self.min_size = min_cluster_size

    def print_new_domains_table(self, table: Dict[str, List[str]],
                                new_domains_table: str) -> None:
        """Prints the new domain `table` in the file `new_domains_table`

        Parameters
        ----------
        table : defaultdict
            table mapping new domain names to list of sequence ids
        new_domains_table : str
            output file for the new domains table
        """
        self.log.info("Printing the new domains table")
        table_file = open(new_domains_table, "w")
        with redirected(stdout=table_file):
            for cluster_name, cluster_list in sorted(table.items()):
                print(cluster_name + "\t" + "\t".join(cluster_list))
        table_file.close()

    def _process_sequence(self, prot: SeqRecord) -> None:
        pass

    def run(self,
            assignment_file: str,
            fasta_file: str,
            interpro_names_file: str,
            clustering_file: str,
            details_file: str = "",
            old_table: str = "",
            prefix: str = "NOVEL") -> None:
        # 1.- set assignment reader and overlap params
        self.assignment_reader = RandomAccessAssignmentReader(assignment_file)

        self.log.info(f"Loading sequences from {fasta_file}...")
        parser = fasta.Parser(open_anything(fasta_file))
        parser = fasta.regexp_remapper(parser,
                                       self.sequence_id_regexp)

        # setup clustering file
        self.clustering = ClusteringFile(self.min_size,
                                         old_table,
                                         prefix=prefix)
        self.clustering_table = self.clustering.process_clustering_file(
            clustering_file)

        # ----------------------------------------

        # loop for each sequence and process it
        for i, seq in enumerate(parser):
            if i % 1000000 == 0:
                self.log.info("Read {} seqs".format(i))
            if len(seq.seq) >= self.maximum:
                self._process_sequence(seq)
