"""Utility classes and routines that do not fit anywhere else"""


class ParseError(Exception):
    """Exception thrown when a parsing error occurred"""

    def __init__(self, msg: str, lineno: int = 1):
        self.parameter = f"{msg} near line {lineno}"

    def __str__(self):
        return repr(self.parameter)
