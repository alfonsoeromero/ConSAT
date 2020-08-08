from gfam.fasta import regexp_remapper
from gfam.fasta.parser import Parser
from gfam.tasks.base import LoggedTask
from gfam.utilities.open_anything import open_anything


class ExtractGeneIDsTask(LoggedTask):
    def process_file(self, filename: str, sequence_id_regexp: str):
        """Processes the given input file"""
        self.log.info("Processing %s..." % filename)

        parser = Parser(open_anything(filename))
        parser = regexp_remapper(parser, sequence_id_regexp)
        for seq in parser:
            print(seq.id)
