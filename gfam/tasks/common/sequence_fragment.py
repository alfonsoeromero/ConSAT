from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SequenceFragment:
    """Represents a fragment of a protein sequence identifier by an id,
        a starting and an ending position"""
    sequence_id: str
    start_pos: int
    end_pos: int

    @classmethod
    def from_str(cls, line: str):
        """Creates a sequence fragment from a string "ID:START-END"
            (same string produced by the __str__ method)

        Parameters
        ----------
        line : str
            String where the `SequenceFragment` is created

        Returns
        -------
        SequenceFragment
            created object
        """
        seq_id, positions = line.split(":")
        start, end = map(int, positions.split("-"))
        return cls(seq_id, start, end)

    def __str__(self) -> str:
        """Gets a string representation of the fragment

        Returns
        -------
        str
            A string representing the fragment: id:start-end
        """
        return f"{self.sequence_id}:{self.start_pos}-{self.end_pos}"

    def __eq__(self, other: object) -> bool:
        """Equal operator, True only if all fields match

        Parameters
        ----------
        other : object
            object to compare with this one, will raise
            `NotImplemented` if it is not a `SequenceFragment`

        Returns
        -------
        bool
            True if both fragments match
        """
        if not isinstance(other, SequenceFragment):
            return NotImplemented
        else:
            return (self.sequence_id == other.sequence_id) and\
                (self.start_pos == other.start_pos) and\
                (self.end_pos == other.end_pos)

    def overlaps(self, other: SequenceFragment) -> bool:
        """Returns if two `SequenceFragment` correspond to the
        same sequence and there is at least one position in common.

        Parameters
        ----------
        other : SequenceFragment
            Sequence fragment to check overlap with.

        Returns
        -------
        bool
            True if this fragment overlaps with `other`.
        """
        return (self.sequence_id == other.sequence_id) and (
            (other.start_pos <= self.start_pos <= other.end_pos) or
            (self.start_pos <= other.start_pos <= self.end_pos))

    def overlap_proportion(self, other: SequenceFragment) -> float:
        """Returns the overlap proportion, i.e.:
        (length of overalapping residues)/(combined length)

        Parameters
        ----------
        other : SequenceFragment
            Sequence fragment to check overlap proportion with.

        Returns
        -------
        float
            a number in [0, 1] measuring the proportion of the overlap.
            If the residues do not overlap the results will be zero.
        """
        if self.overlaps(other):
            min_start = min(self.start_pos, other.start_pos)
            max_end = max(self.end_pos, other.end_pos)
            combined_length = max_end - min_start + 1

            max_start = max(self.start_pos, other.start_pos)
            min_end = min(self.end_pos, other.end_pos)
            common_residues = min_end - max_start + 1

            return common_residues / combined_length
        else:
            return 0.0
