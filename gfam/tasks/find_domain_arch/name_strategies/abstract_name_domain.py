import abc

from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet


class AbstractNameDomainStrategy:
    @abc.abstractmethod
    def get_domain_name(self, fragments: SequenceFragmentSet) -> str:
        """Required Method"""
        pass
