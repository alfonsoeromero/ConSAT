import abc
import re
from collections import defaultdict
from typing import Dict, List, Tuple

from gfam.assignment import SequenceWithAssignments
from gfam.tasks.common.labelled_sequence_fragment_set import \
    LabelledSequenceFragmentSet
from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet
from gfam.tasks.find_domain_arch.name_strategies.abstract_name_domain import \
    AbstractNameDomainStrategy
from gfam.tasks.find_domain_arch.name_strategies.name_domain_with_table import \
    NameDomainWithTable
from gfam.tasks.find_domain_arch.name_strategies.name_domain_without_table import \
    NameDomainWithoutTable


class ClusteringFile:
    """Handles all interactions with the clustering file (output of `cca`) when
        building the architecture in `find_domain_arch`."""

    def __init__(self,
                 min_size: int,
                 old_table: str = "",
                 prefix: str = "NOVEL"):
        """Constructor

        Parameters
        ----------
        min_size : int
            Minimum size (fragments) of a cluster to be considered. Clusters
            with less than this size will be ignored
        old_table : str
            Existing cluster table, if any, to provide names coherent with the previous
            file, by default ''
        prefix : str
            Prefix of the new domain names, by default "NOVEL"
        """

        self.prefix = prefix
        self.min_size = min_size
        self.new_domain_assignments: Dict[str, List] = defaultdict(list)

        self._name_domain_strategy: AbstractNameDomainStrategy
        if old_table:
            self._name_domain_strategy = NameDomainWithTable(
                self.prefix, old_table)
        else:
            self._name_domain_strategy = NameDomainWithoutTable(
                self.prefix)

    def add_new_cluster_assignment_to_sequence(self,
                                               seq: SequenceWithAssignments) -> None:
        """Given a `SequenceWithAssignment` `seq` add all the assignments that
        have new found domains to it, if any

        Parameters
        ----------
        seq : SequenceWithAssignments
            sequence to be completed with new domain assignments
        """
        id = seq.name
        if id in self.new_domain_assignments:
            for start, end, new_domain_id in self.new_domain_assignments[id]:
                seq.assign_(start, end, new_domain_id)

    def process_clustering_file(self, cluster_file: str) -> Dict[str, List[str]]:
        """Given a clustering file

        Parameters
        ----------
        cluster_file : str
            [description]

        Returns
        -------
        [type]
            [description]
        """
        domain_assignment_table: Dict[str, List[str]] = {}

        for line in open(cluster_file):
            fragments = SequenceFragmentSet.from_str(line)
            if fragments.num_different_sequences() < self.min_size:
                continue

            domain_name: str = self._name_domain_strategy.get_domain_name(
                fragments)

            domain_assignment_table[domain_name] = [str(fragment) for fragment in
                                                    fragments]

            for fragment in fragments:
                seq_id = fragment.sequence_id
                self.new_domain_assignments[seq_id].append(
                    (fragment.start_pos, fragment.end_pos, domain_name))

        return domain_assignment_table
