import os
import re
import sqlite3
import tempfile
from re import Pattern


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
        self.regexp: Pattern[str] = re.compile(
            sequence_id_regexp)

        if sequence_id_regexp:
            self._process_protein_id = self._process_protein_id_with_regex
        else:
            self._process_protein_id = self._process_protein_id_without_regex

        self.db_file = self._create_db()
        self.db_handle = self._open_db()

    def _open_db(self) -> sqlite3.Cursor:
        conn = sqlite3.connect(self.db_file)
        return conn.cursor()

    def _create_db(self) -> str:
        db_name = tempfile.NamedTemporaryFile(suffix=".sqlite3", delete=False)
        conn = sqlite3.connect(db_name.name, isolation_level="DEFERRED")
        c = conn.cursor()
        c.execute('''PRAGMA journal_mode = OFF''')
        c.execute("PRAGMA synchronous = OFF")
        c.execute("BEGIN TRANSACTION")

        create_sql = """CREATE TABLE lcr (
           protein_id TEXT,
           left INTEGER,
           right INTEGER
        );"""
        c.execute(create_sql)

        cache = []
        MAX_CACHE_SIZE: int = 100000
        current_prot_id = ""
        for line in open(self.seg_filename):
            if line.startswith(">"):
                current_prot_id = self._process_protein_id(line)
            else:
                fields = line.strip().split()
                left = int(fields[0])
                right = int(fields[2])

                cache.append((current_prot_id, left, right))

                if len(cache) > MAX_CACHE_SIZE:
                    c.executemany(
                        "insert into lcr values (?, ?, ?);", cache)
                    cache = []

        if cache:
            c.executemany("insert into lcr values (?, ?, ?)", cache)

        sql = ("CREATE INDEX index_on_protein ON lcr (protein_id);")
        c.execute(sql)

        # commit the changes to db
        conn.commit()

        # close the connection
        conn.close()
        return db_name.name

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
        self.db_handle.execute(
            "SELECT left, right FROM lcr" +
            f" WHERE protein_id='{protein_id}'")
        return [(int(x[0]), int(x[1])) for x in self.db_handle.fetchall()]

    def delete_temp_file(self) -> None:
        """Delete the associated temporary file"""
        os.unlink(self.db_file)
