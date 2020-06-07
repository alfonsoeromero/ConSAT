import os
import re
import tempfile
from typing import List

import tables

from gfam.tasks.find_unassigned.readers.indexed_assignment_reader import \
    AssignmentPosition


class RandomAccessSEGReader:
    def __init__(self, seg_filename: str,
                 sequence_id_regexp: str = ""):
        """Constructor

        Parameters
        ----------
        seg_filename : str
            full path to the underlying SEG mask file
        sequence_id_regexp : str
            regex to extract the ID from the protein name
            (should contain the <ID> field), default ""
        """
        self.seg_filename = seg_filename
        self.temp_file = tempfile.NamedTemporaryFile(delete=False,
                                                     prefix="seg_reader",
                                                     suffix='.h5')
        if sequence_id_regexp:
            self.regexp = re.compile(sequence_id_regexp)
            self._process_protein_id = self._process_protein_id_with_regex
        else:
            self.regexp = None
            self._process_protein_id = self._process_protein_id_without_regex

        self._build_index()
        self.table = self._open_table()

    def _process_protein_id_without_regex(self, line: str) -> str:
        """Extract the protein ID from a line (get first field and remove >)

        Parameters
        ----------
        line : str
            line that contains the protein id

        Returns
        -------
        str
            extracted protein id
        """
        return line.split()[0][1:]

    def _process_protein_id_with_regex(self, line: str) -> str:
        """Extract the protein ID from a file line using `self.regexp

        Parameters
        ----------
        line : str
            line that contains the protein id

        Returns
        -------
        str
            extracted protein id
        """
        current_prot_id = line.split()[0][1:]
        return self.regexp.sub(r'\g<id>', current_prot_id)

    def _build_index(self) -> None:
        """Builds an inverted index, mapping protein ids
        to list of file positions where the corresponding lines start"""
        inverted_index = tables.open_file(self.temp_file.name, "w")
        grp = inverted_index.create_group("/", "index")
        positions_table = inverted_index.create_table(grp, "index_table",
                                                      AssignmentPosition)

        with open(self.seg_filename, 'r') as in_file:
            prev_offset = -1
            while True:
                line: str = in_file.readline()
                if line:
                    offset: int = in_file.tell()
                    if line.startswith(">"):
                        protein_id = self._process_protein_id(line)
                    else:
                        newrow = positions_table.row

                        newrow["protein_id"] = protein_id
                        newrow["file_position"] = prev_offset
                        newrow.append()
                else:
                    break
                prev_offset = offset

        positions_table.cols.protein_id.create_index()

        inverted_index.flush()
        inverted_index.close()

    def _open_table(self) -> tables.table.Table:
        table = tables.open_file(self.temp_file.name, 'r')
        return table.root.index.index_table

    def _get_postings_list_for_protein(self, protein_id: str) -> List[int]:
        postings = []
        for record in self.table.where(f"protein_id == b'{protein_id}'"):
            postings.append(record["file_position"])
        return postings

    def get_intervals_for_protein(self, protein_id: str):
        """Get the set of intervals (begin, end) from an SEG file
            corresponding to a given protein id.

        Parameters
        ----------
        protein_id : str
            protein identifier that we query for

        Yields
        -------
        Tuple[int, int]
            Interval (begin, end) of the SEG file associated to the
            `protein_id`
        """
        # 1.- get posting list for that particular protein
        postings = self._get_postings_list_for_protein(protein_id)

        # 2.- generator
        with open(self.seg_filename, "r") as in_file:
            for posting in postings:
                in_file.seek(posting)
                line = in_file.readline()
                fields = line.strip().split()
                yield tuple(map(int, [fields[0], fields[2]]))

    def delete_temp_file(self) -> None:
        """Delete the associated temporary file"""
        self.table.close()
        os.unlink(self.temp_file.name)
