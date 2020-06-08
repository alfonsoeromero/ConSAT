import logging
from typing import Dict

from gfam import fasta
from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments
from gfam.interpro import AssignmentReader
from gfam.tasks.base import LoggedTask
from gfam.tasks.find_unassigned.intervals import substract
from gfam.tasks.find_unassigned.low_complexity_regions import \
    read_low_complexity_regions
from gfam.utils import open_anything


class FindUnassignedTask(LoggedTask):
    def __init__(self,
                 max_overlap: int,
                 min_fragment_length: int,
                 min_length: int,
                 sequence_id_regexp: str = None,
                 logger: logging.Logger = None):
        super().__init__(logger)
        AssignmentOverlapChecker.max_overlap = max_overlap
        self.min_fragment_length = min_fragment_length
        self.min_length = min_length
        self.sequence_id_regexp = sequence_id_regexp

        # these are moved here from get_unassigned
        self.regions = []
        self.maximum = None
        # to move somewhere else
        self.seqcat = {}

    def run(self, assignment_file: str, fasta_file: str,
            seg_file: str = None) -> None:

        self.log.info("Processing assignments file...")
        self._process_infile(assignment_file)

        self.log.info("Getting unassigned pieces")
        self._get_unassigned()

        self.log.info("Processing FASTA file...")
        self._process_sequences(fasta_file)

        # if there is a file with low complexity regions, process it
        # and remove low complexity regions from unassigned fragments
        if seg_file is not None and seg_file:
            self.log.info("Removing low complexity regions")
            low_complexity_regions = read_low_complexity_regions(
                seg_file, self.sequence_id_regexp)
            self._remove_low_complexity_regions(low_complexity_regions)

        self.log.info("Getting unassigned pieces")
        self._print_unassigned()

    def _process_sequences(self, fname: str) -> None:
        """Reads and processes a fasta file

        Parameters
        ----------
        fname : str
            path to a fasta file
        """
        self.log.info(f"Loading sequences from {fname}...")
        parser = fasta.Parser(open_anything(fname))
        parser = fasta.regexp_remapper(parser,
                                       self.sequence_id_regexp)
        for i, seq in enumerate(parser):
            if seq.id not in self.seqcat and len(seq.seq) >= self.maximum:
                self.regions.append((seq.id, 1, len(seq.seq)))
            if i % 1000000 == 0:
                self.log.info("Read {} seqs".format(i))
        self.log.info(f"Read a total of {i+1} seqs")

    def _process_infile(self, fname: str, interpro=None) -> None:
        self.log.info("Processing input file: %s" % fname)
        for assignment in AssignmentReader(fname):
            try:
                seq = self.seqcat[assignment.id]
            except KeyError:
                seq = SequenceWithAssignments(assignment.id, assignment.length)
                self.seqcat[assignment.id] = seq
            if seq.length != assignment.length:
                raise ValueError("different lengths encountered "
                                 "for %s: %d and %d" % (seq.name,
                                                        seq.length,
                                                        assignment.length))
            if interpro is not None:
                assignment = assignment.resolve_interpro_ids(interpro)
            seq.assign(assignment, (interpro is not None), interpro=interpro)
        self.log.info(f"Read {len(self.seqcat)} assignements")

    def _get_unassigned(self) -> None:
        self.regions = []
        for seq_id, seq in self.seqcat.items():
            if seq.length < self.min_length:
                continue
            for start, end in seq.unassigned_regions():
                if end-start+1 < self.min_fragment_length:
                    continue
                self.regions.append((seq_id, start, end))
        self.maximum = max(self.min_length,
                           self.min_fragment_length)

    def _remove_low_complexity_regions(self,
                                       low_complexity_regions: Dict) -> None:
        new_regions = []
        for region in self.regions:
            reg_id = region[0]
            if reg_id not in low_complexity_regions:
                new_regions.append(region)
            else:
                lcr_ids = low_complexity_regions[reg_id]
                new_regions.extend(substract(region, lcr_ids))
        self.regions = new_regions
        self.log.info(
            f"{len(self.regions)} remaining regions after LCR removal")

    def _print_unassigned(self):
        maximum = max(self.min_length,
                      self.min_fragment_length)
        for region in self.regions:
            if region[2] - region[1] + 1 >= maximum:
                print("%s\t%d\t%d" % region)
