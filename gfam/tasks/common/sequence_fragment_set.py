from typing import List
from gfam.tasks.common.sequence_fragment import SequenceFragment


class SequenceFragmentSet:
    def __init__(self, fragments: List[SequenceFragment] = []) -> None:
        self.clusters: List[SequenceFragment] = fragments

    @classmethod
    def from_str(cls, line: str) -> 'SequenceFragmentSet':
        fragments = [SequenceFragment.from_str(x)
                     for x in line.strip().split()]
        return cls(fragments)

    def __str__(self) -> str:
        """Returns a string representation of the set

        Returns
        -------
        str
            representation of the set
        """
        return " ".join(map(str, self.clusters))

    def append(self, fragment: SequenceFragment) -> None:
        """Append a `fragment` to the current set

        Parameters
        ----------
        fragment : SequenceFragment
            fragment to append
        """
        self.clusters.append(fragment)

    def size(self) -> int:
        """Return number of fragments of the set

        Returns
        -------
        int
            number of fragments of the set
        """
        return len(self.clusters)

    def num_different_sequences(self) -> int:
        """Number of unique sequence ids in the set

        Returns
        -------
        int
            number of unique sequence ids present in the set
        """
        return len(set(x.sequence_id for x in self.clusters))
