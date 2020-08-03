import sqlite3
import tempfile
from typing import List

from gfam.tasks.find_domain_arch.architecture_assignment import \
    ArchitectureAssignment


class ArchitectureMappingWriter:
    """Responsible for handling a mapping between proteins
        and architectures
    """

    def __init__(self):
        self._BUFFER_SIZE = 100_000
        self._buffer: List[ArchitectureAssignment] = []
        self.db_name = tempfile.NamedTemporaryFile(
            suffix=".sqlite3", delete=False)
        self.conn = sqlite3.connect(
            self.db_name.name, isolation_level="DEFERRED")
        self.c = self.conn.cursor()
        self._create_db()

    def _create_db(self) -> None:
        create_sql = """CREATE TABLE architecture_assignment (
                        protein_id STRING,
                        protein_length INT,
                        architecture STRING,
                        residues_covered INT,
                        architecture_detail STRING,
                        residues_covered_by_novel INT);"""
        self.c.execute(create_sql)
        self.conn.commit()

    def _flush_buffer(self) -> None:
        """Flushes the buffer
        """
        if self._buffer:
            self.c.executemany("""INSERT INTO architecture_assignment
                                  values (?, ?, ?, ?, ?, ?);""",
                               [x.astuple() for x in self._buffer])

    def append(self, assignment: ArchitectureAssignment) -> None:
        self._buffer.append(assignment)
        if len(self._buffer) == self._BUFFER_SIZE:
            self._flush_buffer()
            self._buffer = []

    def finish_insertion(self) -> None:
        """Finishes previous insertion (flush buffer and
        commit insert)"""
        self._flush_buffer()
        self.conn.commit()
        self.conn.close()
