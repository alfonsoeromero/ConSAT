from dataclasses import dataclass

from gfam.tasks.common.sequence_fragment_set import SequenceFragmentSet


@dataclass
class LabelledSequenceFragmentSet:
    """Represent a set of fragments with a label (e.g. a cluster id)"""
    label: str
    fragment_set: SequenceFragmentSet

    @classmethod
    def from_string(cls, line: str) -> 'LabelledSequenceFragmentSet':
        """Build it from a string `line`, where the first token is the
            label and the rest is the `fragment_set`.
        """
        label, rest = line.split(maxsplit=1)
        fragment_set = SequenceFragmentSet.from_str(rest)
        return cls(label, fragment_set)
