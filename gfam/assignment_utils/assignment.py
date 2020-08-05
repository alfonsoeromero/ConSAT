from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Assignment:
    """Class representing a record in an InterPro ``iprscan`` output.

    An InterPro domain assignment has the following fields:

    - ``id``: the ID of the sequence
    - ``length``: the length of the sequence
    - ``start``: the starting position of the region in the sequence
      that is assigned to some InterPro domain (inclusive)
    - ``end``: the ending position (inclusive)
    - ``source``: the assignment source as reported by ``iprscan``
    - ``domain``: the ID of the domain being assigned, according to the
      assignment source
    - ``evalue``: E-value of the assignment if that makes sense,
      ``None`` otherwise
    - ``interpro_id``: the InterPro ID corresponding to ``domain`` in
      ``source``.
    - ``comment``: an arbitrary comment
    """
    id: str
    length: int
    start: int
    end: int
    source: str
    domain: str
    evalue: Optional[float]
    interpro_id: Optional[str]
    comment: Optional[str]

    def get_assigned_length(self):
        """Returns the number of amino acids covered by the assignment
        within the sequence."""
        return self.end - self.start + 1

    def resolve_interpro_ids(self, interpro):
        """If the assignment has an InterPro ID, this method makes sure
        that the domain is equal to the highest common ancestor of the
        InterPro ID in the InterPro tree. If the assignment does not have
        an InterPro ID yet, this method tries to look it up.

        Returns a new tuple which might or might not be equal to this one.
        """
        if self.interpro_id:
            anc = interpro.tree.get_most_remote_ancestor(self.interpro_id)
        else:
            anc = interpro.mapping.get(self.domain)
        if self.domain == anc:
            return self
        return self._replace(domain=anc)

    def short_repr(self):
        """Short representation of this assignment, used in error messages"""
        return "%s(%d-%d)" % (self.domain, self.start, self.end)
