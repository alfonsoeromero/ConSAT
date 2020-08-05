import json
import os
import sqlite3
import tempfile
from dataclasses import asdict
from typing import List

from gfam.interpro import Assignment, AssignmentReader


def _assignment_from_string(s: str) -> Assignment:
    return Assignment(**json.loads(s))


def _assignment_to_string(assignment: Assignment) -> str:
    return json.dumps(asdict(assignment))


class RandomAccessAssignmentReader:
    def __init__(self, filename: str):
        """Access to an assignment file in a random manner

        Parameters
        ----------
        filename : name of the file
            assignment file in InterPro format
        """
        self.filename = filename
        self.db_file = self._create_db()
        self.db_handle = self._open_db()

    def _open_db(self):
        conn = sqlite3.connect(self.db_file)
        return conn.cursor()

    def _create_db(self) -> str:
        db_name = tempfile.NamedTemporaryFile(suffix=".sqlite3", delete=False)
        conn = sqlite3.connect(db_name.name, isolation_level="DEFERRED")
        c = conn.cursor()
        c.execute('''PRAGMA journal_mode = OFF''')
        c.execute("PRAGMA synchronous = OFF")
        c.execute("BEGIN TRANSACTION")

        create_sql = """CREATE TABLE assignments (
           protein_id TEXT,
           assignment TEXT
        );"""
        c.execute(create_sql)

        cache = []
        max_cache_size: int = 100000
        for assignment in AssignmentReader(self.filename):
            cache.append((assignment.id, _assignment_to_string(assignment)))

            if len(cache) > max_cache_size:
                c.executemany("insert into assignments values (?, ?);", cache)
                cache = []

        if cache:
            c.executemany("insert into assignments values (?, ?)", cache)

        sql = ("CREATE INDEX index_on_protein ON assignments (protein_id);")
        c.execute(sql)

        # commit the changes to db
        conn.commit()

        # close the connection
        conn.close()
        return db_name.name

    def assignments_for_protein(self, protein_id: str) -> List[str]:
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
        self.db_handle.execute(
            "SELECT assignment FROM assignments" +
            f" WHERE protein_id='{protein_id}'")
        return [_assignment_from_string(x[0])
                for x in self.db_handle.fetchall()]

    def delete_temp_file(self):
        """Delete the associated temporary file"""
        os.unlink(self.db_file)
