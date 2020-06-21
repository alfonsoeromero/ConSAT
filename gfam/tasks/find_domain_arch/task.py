import logging
from typing import Dict, List

from gfam.tasks.base import LoggedTask
from gfam.utils import redirected


class FindDomainArchTask(LoggedTask):
    def __init__(self,
                 logger: logging.Logger = None):
        super().__init__(logger)

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
