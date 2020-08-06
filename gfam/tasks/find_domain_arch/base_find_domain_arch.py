from gfam.assignment import SequenceWithAssignments
from gfam.config import ConfigurableOptionParser
from gfam.scripts import CommandLineApp
from gfam.scripts.find_unassigned import FindUnassignedApp


class BaseFindDomainArch(CommandLineApp):
    def _add_common_options(self, parser: ConfigurableOptionParser) -> None:
        """Adds the common options of both `find_domain_arch` and
        `find_domain_archi_with_hmms` to the parser

        Parameters
        ----------
        parser : ConfigurableOptionParser
            parser of the application
        """
        parser.add_option("-s", "--min-size", dest="min_size",
                          metavar="N",
                          help="consider only clusters with at least N "
                               "different sequences (not just fragments) "
                               "as novel domains",
                          config_key="analysis:find_domain_arch/"
                                     "min_novel_domain_size",
                          default=2, type=int)
        parser.add_option("-S", "--sequences",
                          dest="sequences_file", metavar="FILE",
                          help="FASTA file containing all the sequences"
                               " of the representative gene model",
                          config_key="file.input.sequences", default=None)
        parser.add_option("-i", "--interpro-parent-child-file",
                          dest="interpro_parent_child_file",
                          metavar="FILE",
                          help="use the InterPro parent-child FILE "
                               "to remap IDs",
                          config_key="file.mapping.interpro_parent_child",
                          default=None)
        parser.add_option("--details",
                          dest="details", metavar="FILE",
                          help="print more details about the domain "
                               "architecture into FILE",
                          config_key="generated/"
                                     "file.domain_architecture_details",
                          default=None)
        parser.add_option("--stats",
                          dest="stats", metavar="FILE",
                          help="print genome-level statistics about "
                               "the domain architectures into FILE",
                          config_key="generated/"
                                     "file.domain_architecture_stats",
                          default=None)
        parser.add_option("-n", "--names",
                          dest="interpro_names_file",
                          metavar="FILE",
                          help="use the given FILE to assign "
                               "InterPro IDs to names",
                          config_key="file.mapping.interpro2name",
                          default=None)
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          help="remap sequence IDs using REGEXP",
                          config_key="sequence_id_regexp",
                          dest="sequence_id_regexp")
        parser.add_option("--max-overlap", metavar="SIZE",
                          help="sets the maximum overlap size allowed between "
                               "assignments of the same data source."
                               " Default: %default",
                          config_key="max_overlap",
                          dest="max_overlap", type=int, default=20)

    def process_interpro_file(self, interpro_file: str,
                              for_hmm_scanning: bool = True) -> None:
        unassigned_app = FindUnassignedApp()
        unassigned_app.set_sequence_id_regexp(self.options.sequence_id_regexp)
        if for_hmm_scanning:
            unassigned_app.process_sequences_file(self.options.sequences_file)
        else:
            unassigned_app.process_sequences_file_old(
                self.options.sequences_file)
        unassigned_app.process_infile(interpro_file, self.interpro)
        self.seqcat = unassigned_app.seqcat
        unassigned_ids = set(unassigned_app.seq_ids_to_length.keys())
        seqcat_ids = set(self.seqcat.keys())
        for seq_id in unassigned_ids - seqcat_ids:
            self.seqcat[seq_id] = SequenceWithAssignments(
                seq_id, unassigned_app.seq_ids_to_length[seq_id])
