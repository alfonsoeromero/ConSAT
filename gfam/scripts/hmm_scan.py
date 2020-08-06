#!/usr/bin/env python
"""HMM scan of a set of pre-trained models in putative protein
domains against a set of sequences
"""

import os
import os.path
import subprocess
import sys
import tempfile
from collections import defaultdict
from operator import itemgetter

from gfam.scripts import CommandLineApp
from gfam.utils import search_file

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2012, Alfonso E. Romero"
__license__ = "GPL"


class HMMScanApp(CommandLineApp):
    """\
    Usage: %prog [options] sequences_file hmm_file

    Given a sequence database in FASTA format, and a HMM file,
    runs hmmsearch an returns the hits of the HMM models in a file
    with the format:

        ID_SEQUENCE:i-j ID_MODEL EVALUE

    where i and j are positions in the sequence, and EVALUE is the
    expected value returned by the HMM tool
    """

    short_name = "hmm_scan"

    def __init__(self, *args, **kwds):
        super(HMMScanApp, self).__init__(*args, **kwds)
        self.hmm_file = None
        self.sequence_file = None

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(HMMScanApp, self).create_parser()

        parser.add_option("-a", dest="num_threads", default=1,
                          type=int, metavar="N",
                          help="instruct hmmsearch to use N cores (threads)."
                               "This option is passed on intact to hmmsearch."
                               " Default: %default",
                          config_key="num_cpu_cores")
        parser.add_option("-e", dest="evalue", default=1e-2,
                          help="evalue threshold for inclusion of sequence"
                               " hits. Default",
                          config_key="analysis:hmm_scan/hmmer_theshold."
                                     " Default: %default")
        parser.add_option("--hmmscan-path", dest="hmmsearch_path",
                          help="uses PATH as the path "
                               "to the hmmsearch executable",
                          config_key="utilities/util.hmmscan",
                          metavar="PATH")
        parser.add_option("--hmmpress-path", dest="hmmpress_path",
                          help="uses PATH as the path "
                               "to the hmmpress executable",
                          config_key="utilities/util.hmmpress",
                          metavar="PATH")
        return parser

    def _create_temp_file(self):
        file_h = tempfile.NamedTemporaryFile(delete=True)
        file_h.close()
        return file_h.name

    def run_real(self):
        """Runs the application and returns the exit code"""
        if not self.args:
            self.args = ["-"]
        else:
            self.sequence_file, self.hmm_file = self.args

        # if the sequences file is empty (i.e. all the sequences are
        # covered by InterPro), then we return
        if os.stat(self.sequence_file).st_size == 0:
            return 0

        # Finds hmmsearch in the current path if needed, ensure that
        # the paths are absolute.
        for util in ["hmmsearch", "hmmpress"]:
            optkey = "%s_path" % util
            path = getattr(self.options, optkey)
            if not path:
                path = search_file(util, executable=True)
                if not path:
                    self.error("cannot find %s in system path" % util)
            setattr(self.options, optkey, os.path.abspath(path))

        # we create a temporal file to store the HMM output
        hmm_output = self._create_temp_file()
        if not self.press_hmm_library():
            return 1
        if not self.run_hmm_search(hmm_output):
            return 1

        table = self.process_hmm_output(hmm_output)
        self.write_table(table)
        return 0

    def write_table(self, table):
        for entry in table:
            print(" ".join(entry))

    def press_hmm_library(self):
        """Runs ``hmmpress`` on the collection of HMM to make
        it a binary file an make it able to be run with hmmsearch
        """
        # First, we check whether the binary files exist
        binary_filenames = [self.hmm_file + ".h3" + i for
                            i in ["f", "i", "m", "p"]]
        if not all([os.path.isfile(fi) for fi in binary_filenames]):
            # the binaries do not exist, we create them
            self.log.info("Pressing the HMM library")
            args = []
            args.extend(["-f"])
            args.extend([self.hmm_file])

            cmd = self.get_tool_cmdline("hmmpress", args)
            if not cmd:
                self.log.fatal("cannot find hmmpress in %s",
                               self.options.hmmpress_path)
                return False

            self.log.info("Executing the following command %s", str(cmd))

            hmmpress_prog = subprocess.Popen(cmd, stdin=open(os.devnull),
                                             stdout=open(os.devnull))
            retcode = hmmpress_prog.wait()

            if retcode != 0:
                self.log.fatal("hmmpress exit code was %d, exiting...",
                               retcode)
                return False

            self.log.info("HMM models pressed")
        return True

    def get_tool_cmdline(self, tool_name, args=None):
        """Given the name of an HMMer tool in `tool_name`, looks up the
        value of the corresponding option from the command line, and tries
        to figure out whether it is a path to a folder containing the tools
        or the path to the tool itself. Returns a list of arguments that
        should be passed to `subprocess.Popen` in order to execute the tool.
        `args` is a set of extra arguments that should be passed to the tool.
         """
        if not args:
            args = []
        tool_path = getattr(self.options, tool_name + "_path")
        if not tool_path:
            tool_path = os.getcwd()

        if os.path.isfile(tool_path):
            # The path exists and points to a real file, so we simply return
            return [tool_path] + args
            # if not os.path.exists(tool_path):
            # The path does not exist. Assume that this is a
            # file referring to the HMMer tools, but the user
            # has the new one. Try to extract the folder and
            # see if the folder exists.
            # _, tool_name = os.path.split(tool_path)

        # Okay, first check for the actual tool in the folder
        full_path = os.path.join(tool_path, tool_name)

        if os.path.isfile(full_path):
            return [full_path] + args
        # Nothing succeeded
        return None

    def process_hmm_line(self, line):
        """Process a line coming from the tabular output of HMMer, taking
        the following fields as a list: sequence, model, evalue, starting
        position of the hit, ending position of the hit, and returning the
        list
        """
        fields = line.split()
        return [fields[0], fields[3], fields[6], fields[19], fields[20]]

    def process_hmm_output(self, hmm_output):
        """Process the tabular output from a HMMer execution, returning a
        list of lists, where each element of the list contains three elements:
        IDsequence:begin-end, HMM_model producing the hit AND score of the hit
        (evalue"""
        hits_per_sequence = defaultdict(list)
        # first, we retrieve all hits per sequence
        for line in open(hmm_output):
            if not line.startswith("#"):
                [seq, model, sco, sfrom, sto] = self.process_hmm_line(line)
                hit = [model, sfrom, sto, sco]
                hits_per_sequence[seq].append(hit)

        # second, an assignment table is built for all sequences,
        # combining hits
        assignment_table = []
        for (sequence, hits) in hits_per_sequence.items():
            sorted_hits = self.sort_hits_by_length(hits)
            accepted_hits = []
            for hit in sorted_hits:
                if not self.overlaps(hit, accepted_hits):
                    accepted_hits.append(hit)
            [real_sequence, limits] = sequence.split(":")
            [base_begin, _] = limits.split("-")

            # here we create the following table
            # ID_SEQUENCE:i-j ID_MODEL EVALUE
            for hit in accepted_hits:
                [hit_model, hit_begin, hit_end, hit_score] = hit
                length = int(hit_end) - int(hit_begin) + 1
                begin = int(base_begin) + int(hit_begin) - 1
                end = begin + length - 1
                new_seq = "{}:{}-{}".format(real_sequence, begin, end)
                assignment_table.append([new_seq,
                                         hit_model,
                                         hit_score])
        return assignment_table

    def overlaps(self, hit, accepted_hits):
        """Checks whether a given hit [model, begin, start, score]
        overlaps with a certain list of hits. The score is not considered
        at all to check the overlap, just the positions
        """
        for accepted_hit in accepted_hits:
            if max(0, min(int(hit[2]), int(accepted_hit[2])) -
                   max(int(hit[1]), int(accepted_hit[1]))) > 0:
                return True
        return False

    def sort_hits_by_length(self, hits):
        """Takes a list of hits [model, begin, start, score]
        and returns it sorted by length
        """
        decorated_hits = [[hit, int(hit[2])-int(hit[1])] for hit in hits]
        decorated_hits.sort(key=itemgetter(1), reverse=True)
        return [hit[0] for hit in decorated_hits]

    def run_hmm_search(self, hmm_output_file):
        """Runs ``hmmsearch`` on the given sequence file, writing the output to the
        given ``hmm_output_file``.
        Returns ``True`` if the execution was successful, ``False`` otherwise.
        """
        self.log.info("Invoking hmmsearch, this might take a long time...")

        args = []
        args.extend(["--cpu", str(self.options.num_threads)])
        args.extend(["--domtblout", hmm_output_file])
        args.extend(["--noali"])
        args.extend(["-o", os.devnull])
        args.extend(["--acc"])
        args.extend(["--incE", str(self.options.evalue)])
        args.append(self.hmm_file)  # hmm database
        args.append(self.sequence_file)  # sequence file

        cmd = self.get_tool_cmdline("hmmsearch", args)
        if not cmd:
            self.log.fatal("cannot find hmmsearch in %s",
                           self.options.hmmscan_path)
            return False

        self.log.info("Executing %s", " ".join(map(str, cmd)))

        hmmscan_prog = subprocess.Popen(cmd, stdin=open(os.devnull),
                                        stdout=open(os.devnull))
        retcode = hmmscan_prog.wait()
        if retcode != 0:
            self.log.fatal("hmmsearch exit code was %d, exiting...", retcode)
            return False

        self.log.info("hmmsearch returned successfully.")
        return True


if __name__ == "__main__":
    sys.exit(HMMScanApp().run())
