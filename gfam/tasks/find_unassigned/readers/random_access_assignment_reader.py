from gfam.interpro import AssignmentReader
from gfam.tasks.find_unassigned.readers.indexed_assignment_reader import \
    IndexedAssignmentReader


class RandomAccessAssignmentReader(AssignmentReader):
    def __init__(self, filename: str):
        """Access to an assignment file in a random manner

        Parameters
        ----------
        filename : name of the file
            assignment file in InterPro format
        """
        super().__init__(filename)
        self.index = IndexedAssignmentReader(filename)

    def assignments_for_protein(self, protein_id: str):
        """Get the list of proteins associated to a given protein
            identifier

        Parameters
        ----------
        protein_id : str
            protein identifier that we query for

        Yields
        -------
        Assignment
            assignments associated to that protein identifier
        """
        for line in self.index.get_lines_with_assignments_for_protein(
                protein_id):
            assignment = self.parse_line(line)
            if assignment is not None:
                yield assignment

    def delete_temp_file(self):
        """Delete the associated temporary file"""
        self.index.delete_temp_file()
