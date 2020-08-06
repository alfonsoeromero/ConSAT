from typing import Optional

from gfam.tasks.assignment_source_filter.exclusion_logger.empty_logger import \
    EmptyLogger
from gfam.tasks.assignment_source_filter.exclusion_logger.file_logger import \
    FileLogger
from gfam.tasks.assignment_source_filter.exclusion_logger.log_strategy import \
    LogStrategy


class ExclusionLogger:
    """Manages the logging of the exclusion reasons of certain sequences"""

    def __init__(self, logger_file: Optional[str] = None):
        self._exclusion_logger: LogStrategy
        if logger_file is None or not logger_file:
            self._exclusion_logger = EmptyLogger()
        else:
            self._exclusion_logger = FileLogger(logger_file)

    def log_exclusion(self, name: str, reason: str) -> None:
        """Adds an entry to the exclusions log file, noting that the
        sequence with the given `name` was excluded from further consideration
        because of the given `reason`. This only works if a `logger_file`
        was provided in the constructor.

        Parameters
        ----------
        name : str
            name of the sequence
        reason : str
            reason for exclusion
        """
        self._exclusion_logger.log_exclusion(name, reason)

    def close(self) -> None:
        """Shuts down logger"""
        self._exclusion_logger.close()
