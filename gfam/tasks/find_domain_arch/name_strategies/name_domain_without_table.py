from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet
from gfam.tasks.find_domain_arch.name_strategies.abstract_name_domain\
    import AbstractNameDomainStrategy


class NameDomainWithoutTable(AbstractNameDomainStrategy):
    """Name new domains in the absence of a domain table"""

    def __init__(self, prefix: str):
        """Constructor

        Parameters
        ----------
        prefix : str
            Prefix for new domains
        """
        self.current_cluster_id = 0
        self.prefix = prefix

    def get_domain_name(self, fragments: SequenceFragmentSet) -> str:
        """Gets the domain name

        Parameters
        ----------
        fragments : SequenceFragmentSet
            fragment set that will compose the new domain (unused in this
            case since there is no previous table)

        Returns
        -------
        str
            name of the new domain, usually prefix followed by 5 digits
            starting from 0
        """
        cluster_id: int = self.current_cluster_id
        domain_name = f"{self.prefix}{cluster_id:05}"
        self.current_cluster_id += 1
        return domain_name
