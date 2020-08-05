from dataclasses import dataclass


@dataclass
class Slice:
    """Represents a slice of a given sequence between two positions"""
    sequence_id: str
    left: int
    right: int
