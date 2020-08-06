from enum import Enum


class OverlapType(Enum):
    """Enum describing the different overlap types that can be detected
    by `AssignmentOverlapChecker`. See the documentation of
    `AssignmentOverlapChecker.check_single` for more details.
    """
    NO_OVERLAP: str = "NO_OVERLAP"
    INSERTION: str = "INSERTION"
    INSERTION_DIFFERENT: str = "INSERTION_DIFFERENT"
    SYNONYM_INSERTION: str = "SYNONYM_INSERTION"
    DIFFERENT: str = "DIFFERENT"
    OVERLAP: str = "OVERLAP"
