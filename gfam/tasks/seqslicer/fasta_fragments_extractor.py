import logging
import sys
from typing import Dict, List, Optional, Tuple

from gfam import fasta
from gfam.sequence import SeqRecord
from gfam.tasks.seqslicer.slice import Slice
from gfam.utilities.open_anything import open_anything


class FastaFragmentsExtractor:
    def __init__(self, output_file: Optional[str] = None,
                 sequence_id_regexp: Optional[str] = None,
                 try_alternative_splicing: bool = False,
                 keep_ids: bool = False,
                 log: logging.Logger = None):
        self.output_file = output_file
        self.sequence_id_regexp = sequence_id_regexp
        self.try_alternative_splicing = try_alternative_splicing
        self.keep_ids = keep_ids
        if log is None:
            self.log = logging.getLogger(__name__)
        else:
            self.log = log

    def _setup_writers(self) -> List:
        writer = fasta.FastWriter(sys.stdout)
        writers = [writer]
        if self.output_file is not None:
            output_fd = open(self.output_file, "w")
            writer_file = fasta.FastWriter(output_fd)
            writers.append(writer_file)
        return writers

    def _get_seq_id_from_sequence(self, seq: SeqRecord,
                                  dict_slices: Dict[str, List[Slice]]) ->\
            Tuple[str, bool]:
        seq_id = seq.id
        is_found = True
        if seq_id not in dict_slices:
            if self.try_alternative_splicing:
                seq_id = seq_id.strip().rstrip(".1")
                if seq_id not in dict_slices:
                    is_found = False
            else:
                is_found = False
        return seq_id, is_found

    def process_sequences_file(self, seq_file: str,
                               dict_slices: Dict[str, List[Slice]]) -> int:
        """Processes the sequences one by one, extracting all the pieces into
        an output fasta file"""

        parser = fasta.Parser(open_anything(seq_file))
        parser = fasta.regexp_remapper(parser, self.sequence_id_regexp)

        ids_to_process = set(dict_slices.keys())

        writers: List = self._setup_writers()

        for seq in parser:
            seq_id, is_found = self._get_seq_id_from_sequence(seq, dict_slices)
            if not is_found:
                continue

            sequence = seq.seq
            length_seq = len(sequence)
            ids_to_process.remove(seq_id)

            for seq_slice in dict_slices[seq_id]:
                left, right = seq_slice.left, seq_slice.right

                right = min(right, length_seq)
                # just in case...

                if left > right:
                    # again, just in case
                    self.log.warning("Problem with fragment of %s, "
                                     "the right part is smaller than "
                                     "the left", seq_id)
                    continue

                new_record = None

                if left == 1 and right == length_seq:
                    new_record = seq.fragment(not self.keep_ids)
                else:
                    if not self.keep_ids:
                        new_id = "%s:%d-%d" % (seq_id, left, right)
                    else:
                        new_id = seq_id
                    new_record = SeqRecord(sequence[(left-1):right],
                                           id=new_id,
                                           name=seq.name,
                                           description="")
                for writer in writers:
                    writer.write(new_record)

        if len(writers) > 1:
            writers[1].close()

        if ids_to_process:
            self.log.critical("The following identifiers of sequences (%s)"
                              " were found in the fragments file, but not in "
                              " the fasta file ", ",".join(ids_to_process))
            return 1
        return 0
