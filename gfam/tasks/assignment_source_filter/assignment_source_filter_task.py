import logging
from collections import defaultdict
from typing import Set, Union

from gfam.tasks.assignment_source_filter.assignment_filter import \
    AssignmentFilter
from gfam.tasks.assignment_source_filter.assignment_reader_with_filters import\
    AssignmentReaderWithFilters
from gfam.utilities.complementer_set import ComplementerSet


class AssignmentSourceFilterTask:
    def __init__(self,
                 assignment_reader_with_filters: AssignmentReaderWithFilters,
                 assignment_filter: AssignmentFilter,
                 log: logging.Logger):
        self.assignment_reader = assignment_reader_with_filters
        self.assignment_filter = assignment_filter
        self.log = log

    def process_infile(self,
                       valid_sequence_ids: Union[ComplementerSet, Set[str]]):
        self.log.info("Processing assignments...")

        current_id = None
        assignments_by_source: defaultdict = defaultdict(list)
        for assignment, line in self.assignment_reader.assignments_and_lines():
            if assignment.id != current_id:
                if current_id is not None and\
                        current_id in valid_sequence_ids:
                    self.assignment_filter.filter_and_print_assignments(
                        current_id, assignments_by_source)
                current_id = assignment.id
                assignments_by_source = defaultdict(list)

            assignments_by_source[assignment.source].append((assignment, line))

        # ...and the last batch
        if current_id is not None and current_id in valid_sequence_ids:
            self.assignment_filter.filter_and_print_assignments(
                current_id, assignments_by_source)
