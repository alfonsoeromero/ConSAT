from gfam.tasks.assignment_source_filter.exclusion_logger.log_strategy import\
    LogStrategy


class FileLogger(LogStrategy):
    """Logs messages to a file"""

    def __init__(self, logger_file: str):
        self._log_file = open(logger_file, "a+")

    def log_exclusion(self, name: str, reason: str) -> None:
        self._log_file.write(f"{name}: {reason}\n")

    def close(self) -> None:
        self._log_file.close()
