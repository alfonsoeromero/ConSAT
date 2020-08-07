import logging
from typing import Optional, Set, Union

from gfam.utilities.complementer_set import ComplementerSet
from gfam.utilities.open_anything import open_anything


class ValidSequenceIdsFactory:
    """Reads the list of valid ids from a file, if specified
    """

    def __init__(self, log: logging.Logger, file_ids: Optional[str] = None):
        self.log = log
        self.file_ids = file_ids

    def get(self) -> Union[ComplementerSet, Set]:
        """Get list of valid ids (or nothing if the file was not given)

        Returns
        -------
        Union[ComplementerSet, Set]
            list of valid ids
        """
        if self.file_ids is None:
            return ComplementerSet()
        else:
            self.log.info(f"Loading sequence IDs from {self.file_ids}...")
            return set([line.strip()
                        for line in open_anything(self.file_ids)])
