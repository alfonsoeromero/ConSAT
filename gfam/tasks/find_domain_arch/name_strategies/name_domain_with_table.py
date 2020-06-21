import re
from collections import defaultdict
from typing import Dict

from gfam.tasks.common.labelled_sequence_fragment_set import \
    LabelledSequenceFragmentSet
from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet
from gfam.tasks.find_domain_arch.domain_name_strategies.abstract_name_domain import \
    AbstractNameDomainStrategy


class NameDomainWithTable(AbstractNameDomainStrategy):
    """Name new domains in the presence of an existing domain table"""

    def __init__(self, prefix: str, table: str):
        self.prefix = prefix
        self.current_cluster_id = -1
        self._process_existing_cluster_table(table)

    def _find_domain_id(self, fragments: SequenceFragmentSet) -> str:
        """Maps a set of fragments to the most likely
        cluster from the old_cluster_trable, if this is possible.
        The identifier (number) of this cluster is returned, if
        anyone is found, or -1 if no matching cluster is found.
        """
        # 1.- we vote each possible cluster
        votes_per_cluster: Dict[str, int] = defaultdict(int)
        for fragment in fragments:
            s_fragment = str(fragment)
            if s_fragment in self.cluster_per_fragment:
                # 1.1.- if the fragment has a previously assigned cluster,
                # we vote on it
                voted_cluster: str = self.cluster_per_fragment[s_fragment]
                votes_per_cluster[voted_cluster] += 1
            else:
                # 1.2.- if not, we check for fragments near this one in the
                # same sequence
                sequence = fragment.sequence_id

                if sequence in self.fragments_per_seq:
                    # If there are other fragments in the same sequence
                    # we search for the one with the highes overlap...
                    match_scores: Dict[str, float] = {
                        str(other): fragment.overlap_proportion(other)
                        for other in self.fragments_per_seq[sequence]
                        if fragment.overlaps(other)
                    }

                    if match_scores:
                        matched = max(match_scores,
                                      key=match_scores.get)
                        cluster_id = self.cluster_per_fragment[matched]
                        votes_per_cluster[cluster_id] += 1

        # 2.- we count the votes and take a decision
        if votes_per_cluster:
            return max(votes_per_cluster, key=votes_per_cluster.get)
        else:
            return ""

    def _process_existing_cluster_table(self, table_file: str):
        """Process an existing `table_file`, i.e., a file with one
        line per cluster where the first line token is the cluster id
        and the rest of the tokens are the sequence fragments.

        It builds 3 structures: `self.fragments_per_cluster`, the set
        of fragments assigned to each cluster, `self.cluster_per_fragment`
        (the inverse mapping), and `self.fragments_per_seq`, the list of
        fragments for a given sequence.

        Parameters
        ----------
        table_file : str
            previous table file
        """
        self.fragments_per_cluster: defaultdict = defaultdict(list)
        self.cluster_per_fragment: Dict[str, str] = {}
        self.fragments_per_seq = defaultdict(list)

        for line in open(table_file):
            cluster_set = LabelledSequenceFragmentSet.from_string(line)
            cluster_name = cluster_set.label
            fragments = cluster_set.fragment_set
            cluster_num = int(re.findall(r'\d+', cluster_name)[0])
            self.current_cluster_id = max(self.current_cluster_id, cluster_num)

            self.fragments_per_cluster[cluster_name] = fragments
            for fragment in fragments:
                sequence = fragment.sequence_id
                self.cluster_per_fragment[str(fragment)] = cluster_name
                self.fragments_per_seq[sequence].append(fragment)

    def get_domain_name(self, fragments: SequenceFragmentSet) -> str:
        domain_id = self._find_domain_id(fragments)
        if domain_id:
            domain_name = domain_id
        else:
            domain_name = self.prefix +\
                "%05d" % self.current_cluster_id
            self.current_cluster_id += 1
        return domain_name
