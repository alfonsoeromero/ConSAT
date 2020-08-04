import logging
from collections import defaultdict
from typing import Set, Union

from gfam.assignment import EValueFilter, SequenceWithAssignments
from gfam.interpro import AssignmentReader
from gfam.tasks.assignment_source_filter.exclusion_logger.exclusion_logger\
    import ExclusionLogger
from gfam.utils import complementerset


class AssignmentSourceFilterTask:
    def __init__(self, ignored: Set[str],
                 exclusion_logger: ExclusionLogger,
                 valid_sequence_ids: Union[complementerset, Set[str]],
                 max_e: float,
                 stages_from_config,
                 interpro,
                 log: logging.Logger):
        self.ignored = ignored
        self.log = log
        self.max_e = max_e
        self.exclusion_logger = exclusion_logger
        self.stages_from_config = stages_from_config
        self.valid_sequence_ids = valid_sequence_ids
        self.interpro = interpro

    def process_infile(self, fname: str):
        self.log.info(f"Processing {fname}...")

        current_id = None
        assignments_by_source: defaultdict = defaultdict(list)
        evalue_filter = EValueFilter.FromString(self.max_e)

        reader = AssignmentReader(fname)
        for assignment, line in reader.assignments_and_lines():
            if assignment.id != current_id:
                self.filter_and_print_assignments(current_id,
                                                  assignments_by_source)
                current_id = assignment.id
                assignments_by_source = defaultdict(list)

            if (assignment.source in self.ignored) or\
                (assignment.evalue is not None and
                 not evalue_filter.is_acceptable(assignment)):
                continue
            assignments_by_source[assignment.source].append((assignment, line))

        # ...and the last batch
        self.filter_and_print_assignments(current_id, assignments_by_source)

    def filter_assignments(self, name, assignments_by_source):
        """Given a sequence name and its assignments ordered in a dict by
        their sources, selects a representative assignment set based on the
        rules outlined in the documentation of `FindUnassignedApp`.
        """

        if not assignments_by_source:
            self.exclusion_logger.log_exclusion(
                name, "no assignments in the input data file " +
                "passed the filters")
            return []

        # Determine the length of the sequence (and check that the length is
        # the same across all assignments; if not, then the input file is
        # inconsistent and the sequence will be skipped).
        source = list(assignments_by_source.keys())[0]
        seq_length = assignments_by_source[source][0][0].length
        for _source, assignments in assignments_by_source.items():
            if any(assignment.length != seq_length
                   for assignment, _ in assignments):
                self.log.warning("Sequence %s has multiple assignments with "
                                 "different sequence lengths in the "
                                 "input file, skipping" % name)
                self.exclusion_logger.log_exclusion(
                    name, "ambiguous sequence length in input file")
                return []

        # Initially, the result is empty
        result = []

        # Set up the stages
        stages = self.stages_from_config

        # The first stage is treated specially as we have to select a single
        # source thas has the largest coverage. In the remaining stages, we
        # are allowed to cherrypick from different sources.

        # First, find the data source which covers the most of the sequence
        # and is allowed in stage 1
        first_stage = stages[0]
        coverage = {}
        for source, assignments in assignments_by_source.items():
            # Exclude those sources that we don't consider in the first stage
            if source not in first_stage:
                continue

            # Calculate the coverage: we add all the residues covered by
            # each sequence, not taking overlaps into consideration (by the
            # moment)
            seq = SequenceWithAssignments(name, seq_length)
            for a, _ in assignments:
                seq.assign(a, False, interpro=self.interpro)
            coverage[source] = seq.coverage()

        # Find the source giving the best coverage, add its domains into
        # the current assignment.
        seq = SequenceWithAssignments(name, seq_length)
        if coverage:
            best_source = max(coverage.keys(), key=coverage.__getitem__)
            sorted_assignments = sorted(assignments_by_source[best_source],
                                        key=lambda x:
                                            x[0].get_assigned_length(),
                                        reverse=True)
            for a, line in sorted_assignments:
                line = line.strip()
                if seq.assign(a, True, interpro=self.interpro):
                    tab_count = list(line).count("\t")
                    if tab_count < 13:
                        line = line + "\t" * (13-tab_count)
                    result.append("%s\t%s" % (line, 1))
        else:
            best_source = None

        # Collect the unused assignments (not from the best source)
        # into unused_assignments
        unused_assignments = []
        for source, assignments in assignments_by_source.items():
            if source == best_source:
                continue
            unused_assignments.extend(assignments)

        if not unused_assignments:
            return result

        # Try to fill the unassigned regions with the rest of the assignments
        # that were unused so far, starting from the longest assignment.
        unused_assignments.sort(key=lambda x: -x[0].get_assigned_length())

        # Okay, we're done with the first stage, process the rest.
        # idx_to_stage will contain the indices of the selected
        # assignments as keys and the number of the corresponding
        # stage in which they were selected as values.
        idx_to_stage = {}
        for stage_no, sources in enumerate(stages[1:]):
            for idx, (a, _) in enumerate(unused_assignments):
                if a.source in sources and seq.assign(a, True,
                                                      interpro=self.interpro):
                    idx_to_stage[idx] = stage_no+2
        for idx in sorted(idx_to_stage.keys()):
            row = unused_assignments[idx][1].strip()
            tab_count = list(row).count("\t")
            if tab_count < 13:
                row = row + "\t" * (13-tab_count)
            result.append("%s\t%s" % (row, idx_to_stage[idx]))

        if not result:
            self.exclusion_logger.log_exclusion(
                name, "no assignments were selected after "
                "executing all the stages")

        return result

    def filter_and_print_assignments(self, name, assignments_by_source):
        """Filters and prints the list of assignments of the gene with the
        given `name`. `assignments_by_source` must contain the list of
        domain assignments, sorted by data source."""
        if name is None:
            return
        if name not in self.valid_sequence_ids:
            self.exclusion_logger.log_exclusion(
                name, "not in the list of valid gene IDs")
            return
        for row in self.filter_assignments(name, assignments_by_source):
            print(row)
