from gfam.tasks.assignment_source_filter.exclusion_logger.log_strategy import\
    LogStrategy


class EmptyLogger(LogStrategy):
    """Empty logger, does nothing, used when no
        logger file is provided"""

    def __init__(self):
        pass

    def log_exclusion(self, name: str, log_message: str) -> None:
        """Does nothing"""
        pass
