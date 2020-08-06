import logging
from typing import List, DefaultDict, Optional

from gfam.assignment import SequenceWithAssignments
from gfam.interpro import InterPro
from gfam.tasks.assignment_source_filter.exclusion_logger.exclusion_logger\
    import ExclusionLogger


class AssignmentFilter:
    def __init__(self, exclusion_logger: ExclusionLogger,
                 stages_from_config,
                 interpro: InterPro,
                 log: logging.Logger):
        self.exclusion_logger = exclusion_logger
        self.stages_from_config = stages_from_config
        self.interpro = interpro
        self.log = log
        self.seq_length: int = 0

    def filter_and_print_assignments(self, name: str,
                                     assignments_by_source:
                                         DefaultDict[str, List]) -> None:
        """Filters and prints the list of assignments of the gene with the
            given `name`. `assignments_by_source` must contain the list of
            domain assignments, sorted by data source."""

        if not assignments_by_source:
            self.exclusion_logger.log_exclusion(
                name, "no assignments in the input data file " +
                "passed the filters")
        else:
            results = self._filter_assignments(name, assignments_by_source)
            if not results:
                self.exclusion_logger.log_exclusion(
                    name, "no assignments were selected after "
                    "executing all the stages")
            else:
                for row in results:
                    print(row)

    def _assignments_are_of_different_length(self,
                                             assignments_by_source,
                                             name) -> bool:
        # Determine the length of the sequence (and check that the length is
        # the same across all assignments; if not, then the input file is
        # inconsistent and the sequence will be skipped).

        source = list(assignments_by_source.keys())[0]
        self.seq_length = assignments_by_source[source][0][0].length
        for _source, assignments in assignments_by_source.items():
            if any(assignment.length != self.seq_length
                   for assignment, _ in assignments):
                self.log.warning(f"Sequence {name} has multiple assignments "
                                 "with different sequence lengths in the "
                                 "input file, skipping")
                self.exclusion_logger.log_exclusion(
                    name, "ambiguous sequence length in input file")
                return True
        return False

    def _get_coverage_for_assignments(self, name: str,
                                      assignments) -> int:
        seq = SequenceWithAssignments(name, self.seq_length)
        for a, _ in assignments:
            seq.assign(a, False, interpro=self.interpro)
        return seq.coverage()

    def _complete_line_until_n_tabs(self, line: str, n_tabs: int = 13) -> str:
        tab_count = list(line).count("\t")
        if tab_count < n_tabs:
            return line + "\t" * (n_tabs-tab_count)
        else:
            return line

    def _filter_assignments(self, name: str,
                            assignments_by_source: DefaultDict[str, List])\
            -> List[str]:
        """Given a sequence name and its assignments ordered in a dict by
        their sources, selects a representative assignment set based on the
        rules outlined in the documentation of `FindUnassignedApp`.
        """
        if self._assignments_are_of_different_length(assignments_by_source,
                                                     name):
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
            coverage[source] = self._get_coverage_for_assignments(
                name, assignments)

        # Find the source giving the best coverage, add its domains into
        # the current assignment.
        seq = SequenceWithAssignments(name, self.seq_length)
        best_source: Optional[str]
        if coverage:
            best_source = max(coverage.keys(), key=coverage.__getitem__)
            sorted_assignments = sorted(assignments_by_source[best_source],
                                        key=lambda x:
                                            x[0].get_assigned_length(),
                                        reverse=True)
            for a, line in sorted_assignments:
                if seq.assign(a, True, interpro=self.interpro):
                    line = self._complete_line_until_n_tabs(line.strip())
                    result.append(f"{line}\t1")
        else:
            best_source = None

        # Collect the unused assignments (not from the best source)
        # into unused_assignments
        unused_assignments = []
        for source, assignments in assignments_by_source.items():
            if source != best_source:
                unused_assignments.extend(assignments)

        if not unused_assignments:
            return result

        # Try to fill the unassigned regions with the rest of the assignments
        # that were unused so far, starting from the longest assignment.
        unused_assignments.sort(
            key=lambda x: x[0].get_assigned_length(), reverse=True)

        self._fill_with_unused_assignments(unused_assignments,
                                           result, seq, stages)
        return result

    def _fill_with_unused_assignments(self, unused_assignments,
                                      result, seq, stages) -> None:
        """Try to fill the unassigned regions with the rest of the assignments
        that were unused so far, starting from the longest assignment."""
        unused_assignments.sort(
            key=lambda x: x[0].get_assigned_length(), reverse=True)

        # Okay, we're done with the first stage, process the rest.
        # idx_to_stage will contain the indices of the selected
        # assignments as keys and the number of the corresponding
        # stage in which they were selected as values.
        idx_to_stage = {}
        for stage_no, sources in enumerate(stages[1:]):
            for idx, (a, _) in enumerate(unused_assignments):
                if a.source in sources and seq.assign(a, True,
                                                      interpro=self.interpro):
                    idx_to_stage[idx] = stage_no + 2

        for idx in sorted(idx_to_stage.keys()):
            row = self._complete_line_until_n_tabs(
                unused_assignments[idx][1].strip())
            result.append(f"{row}\t{idx_to_stage[idx]}")
