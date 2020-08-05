import array
import logging
from collections import defaultdict
from typing import Optional

from gfam.modula.log import get_logger
from gfam.tasks.seqslicer.slice_builder import SliceBuilder
from gfam.utils import open_anything


class SliceFile:
    def get_slice_file_as_dict(self, slice_file: str) -> defaultdict:
        """Gets the slice file as a dictionary of lists"""
        slice_dict = defaultdict(list)

        for line in open_anything(slice_file):
            if not line.strip():
                continue
            sequence_slice = SliceBuilder.from_string(line)
            slice_dict[sequence_slice.sequence_id].append(sequence_slice)
        return slice_dict
