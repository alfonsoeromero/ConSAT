import logging

from gfam import fasta
from gfam.assignment import AssignmentOverlapChecker, SequenceWithAssignments

from gfam.interpro import Assignment
from gfam.sequence import SeqRecord
from gfam.tasks.find_unassigned.readers.random_access_assignment_reader import\
    RandomAccessAssignmentReader
from gfam.tasks.find_unassigned.readers.random_access_seg_reader import \
    RandomAccessSEGReader
from gfam.tasks.base import LoggedTask
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

        self.maximum = max(self.min_length,
                           self.min_fragment_length)

    def _preload_readers(self, assignment_file: str, seg_file: str) -> None:
        self.assignment_reader = RandomAccessAssignmentReader(assignment_file)
        if seg_file is not None and seg_file:
            self.seg_reader = RandomAccessSEGReader(seg_file,
                                                    self.sequence_id_regexp)

    def _remove_temp_files(self) -> None:
        self.assignment_reader.delete_temp_file()
        if self.seg_reader:
            self.seg_reader.delete_temp_file()

    def _process_sequence(self, prot: SeqRecord) -> None:
        protein_id = prot.id
        protein_length = len(prot.seq)

        # 0.- create a sequence with assignments
        seq = SequenceWithAssignments(protein_id,
                                      protein_length)

        # 1.- get assignments and low complexity regions for the protein
        assignments_for_sequence = list(self.assignment_reader.
                                        assignments_for_protein(protein_id))
        lcrs_for_sequence = [] if self.seg_reader is None else\
            self.seg_reader.get_intervals_for_protein(protein_id)

        # 2.- assign every assignment and LCR interval, if any
        for assignment in assignments_for_sequence:
            seq.assign(assignment, False)

        for interval in lcrs_for_sequence:
            interval_assignment = Assignment(protein_id, protein_length,
                                             interval[0], interval[1], "SEG",
                                             "SEG", 0.0, None, "")
            seq.assign(interval_assignment, False)

        # 3.- print unassigned regions
        for start, end in seq.unassigned_regions():
            region_length = end - start + 1
            if region_length >= self.maximum:
                print(f"{protein_id}\t{start}\t{end}")

    def run(self, assignment_file: str, fasta_file: str,
            seg_file: str = "") -> None:
        # preloads the reader
        self._preload_readers(assignment_file, seg_file)

        self.log.info(f"Loading sequences from {fasta_file}...")
        parser = fasta.Parser(open_anything(fasta_file))
        parser = fasta.regexp_remapper(parser,
                                       self.sequence_id_regexp)
        # loop for each sequence and process it
        for i, seq in enumerate(parser):
            if i % 1000000 == 0:
                self.log.info("Read {} seqs".format(i))
            if len(seq.seq) >= self.maximum:
                self._process_sequence(seq)

        self.log.info(f"Read a total of {i+1} seqs")

        # remove files
        self._remove_temp_files()
