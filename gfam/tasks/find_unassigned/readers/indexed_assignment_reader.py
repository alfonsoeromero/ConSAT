import os
import tempfile
from typing import List

import tables


class AssignmentPosition(tables.IsDescription):
    """Pytables registry containing protein_id and assignment
    file offset"""
    protein_id = tables.StringCol(itemsize=255, pos=0)
    file_position = tables.Int64Col(pos=1)


class IndexedAssignmentReader:
    def __init__(self, interpro_filename: str):
        """Constructor

        Parameters
        ----------
        interpro_filename : str
            full path to the underlying InterPro file
        """
        self.interpro_filename = interpro_filename
        self.temp_file = tempfile.NamedTemporaryFile(
            delete=False, prefix="assignment_reader",
            suffix='.h5')
        self._build_index()
        self.table = self._open_table()

    def _get_num_lines_file(self, fname: str) -> int:
        return sum(1 for line in open(fname))

    def _build_index(self) -> None:
        """Builds an inverted index, mapping protein ids
        to list of file positions where the corresponding lines start"""
        inverted_index = tables.open_file(self.temp_file.name, "w")
        grp = inverted_index.create_group("/", "index")
        num_expected_rows = self._get_num_lines_file(self.interpro_filename)

        positions_table = inverted_index.create_table(
            grp, "index_table", AssignmentPosition,
            expectedrows=num_expected_rows)

        with open(self.interpro_filename, 'r') as in_file:
            while True:
                offset: int = in_file.tell()
                line: str = in_file.readline()
                if not line:
                    break
                protein_id = line.split("\t", 1)[0]
                newrow = positions_table.row

                newrow["protein_id"] = protein_id
                newrow["file_position"] = offset
                newrow.append()

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

    def get_lines_with_assignments_for_protein(self, protein_id: str):
        """Get the set of lines from an InterPro file corresponding to a
            given protein id.

        Parameters
        ----------
        protein_id : str
            protein identifier that we query for

        Yields
        -------
        str
            Lines of the InterPro file associated to the `protein_id`
        """
        # 1.- get posting list for that particular protein
        postings = self._get_postings_list_for_protein(protein_id)

        # 2.- generator
        with open(self.interpro_filename, "r") as in_file:
            for posting in postings:
                in_file.seek(posting)
                yield in_file.readline()

    def delete_temp_file(self) -> None:
        """Delete the associated temporary file"""
        os.unlink(self.temp_file.name)
