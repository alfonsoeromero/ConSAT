from abc import ABC, abstractmethod


class LogStrategy(ABC):
    """Abstract logging strategy"""
    @abstractmethod
    def log_exclusion(self, name: str, reason: str) -> None:
        """Logs the reason message to exclude a given sequence

        Parameters
        ----------
        name : str
            name of the excluded sequence

        reason : str
            reason for which sequence `name` was excluded
        """

    def close(self) -> None:
        """Shuts down logger"""
        pass
