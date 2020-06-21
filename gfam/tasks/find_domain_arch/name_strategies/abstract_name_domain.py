import abc

from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet


class AbstractNameDomainStrategy:
    def __init__(self):
        self.current_cluster_id = 0

    @abc.abstractmethod
    def get_domain_name(self, fragments: SequenceFragmentSet) -> str:
        """Required Method"""
        pass
