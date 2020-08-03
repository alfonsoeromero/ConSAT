import abc
import sqlite3


class AbstractArchitectureReport(abc.ABCMeta):
    def __init__(self, db_name: str,
                 buffer_limit: int = 100_000):
        """Constructor

        Parameters
        ----------
        db_name : str
            Database name (i.e. sqlite3 file name)
        buffer_limit : int, optional
            size of the internal buffer, that is, max. no. of
            rows to keep in memory from the table, by default 100_000.
        """
        self.conn = sqlite3.connect(
            self.db_name.name, isolation_level="DEFERRED")
        self.limit: int = buffer_limit

    @abc.abstractmethod
    def _get_query(self) -> str:
        """Gets the query corresponding to a particular report
        (abstract method).

        Returns
        -------
        str
            SELECT query returned.
        """
        pass

    def get_report(self, output_file: str, sep: str = "\t") -> None:
        """Gets a particular report as a `sep`-separated file
        called `output_file`.

        Parameters
        ----------
        output_file : str
            file to write the report at.
        sep: str
            separator character between field, default = '\t'
        """
        query = self._get_query()

        c = self.conn.execute(query)

        with open(output_file, "w") as out:
            while c:
                result = c.fetchmany(size=self.buffer_limit)
                for res in result:
                    out.write(f"{sep.join(res)}\n")
