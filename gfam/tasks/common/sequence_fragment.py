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
