#!/usr/bin/env python
"""All-against-all BLASTing of a given set of sequences.

"""

from __future__ import with_statement

import os
import os.path
import subprocess
import sys
import tempfile

from gfam import fasta
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything, search_file, temporary_dir

__author__  = "Alfonso E. Romero"
__email__   = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2012, Alfonso E. Romero"
__license__ = "GPL"

class HMMScanApp(CommandLineApp):
    """\
    Usage: %prog [options] sequences_file hmm_file

    Given a sequence database in FASTA format, and a HMM file, 
    runs hmm_scan an returns the hits of the HMM models in a file
    with the format:

    	ID_SEQUENCEi:j ID_MODEL EVALUE

   where i and j are positions in the sequence, and EVALUE is the
   expected value returned by the HMM tool    
    """

    short_name = "hmm_scan"

    def create_parser(self):
        """Creates the command line parser"""
        parser = super(HMMScanApp, self).create_parser()

        parser.add_option("-a", dest="num_threads", default=1,
                type=int, metavar="N",
                help="instruct hmm_scan to use N cores. This option "
                "is passed on intact to hmm_scan. Default: %default",
                config_key="num_cpu_cores"
        )
        parser.add_option("-e", dest="evalue", default=1e-2,
                help="evalue threshold for inclusion of sequence hits. Default",
		config_key="analysis:hmm_scan/hmmer_theshold. Default: %default"
        )
        parser.add_option("--hmmscan-path", dest="hmmscan_path",
                help="uses PATH as the path to the hmmscan executable",
                config_key="utilities/util.hmmscan", metavar="PATH"
        )
        parser.add_option("--hmmpress-path", dest="hmmpress_path",
                help="uses PATH as the path to the hmmpress executable",
                config_key="utilities/util.hmmpress", metavar="PATH"
        )
        return parser

    def create_temp_file(self):
        f = tempfile.NamedTemporaryFile(delete=True)
	print "Temporal file ", f.name, " created"
	f.close()
	return f.name

    def run_real(self):
        """Runs the application and returns the exit code"""
        if not self.args:
            self.args = ["-"]
        else:
	    self.sequence_file, self.hmm_file = self.args

        # Find hmm_scan in the current path if needed, ensure that
        # the paths are absolute.
        for util in ["hmmscan", "hmmpress"]:
            optkey = "%s_path" % util
            path = getattr(self.options, optkey)
            if not path:
                path = search_file(util, executable=True)
                if not path:
                    self.error("cannot find %s in system path" % util)
            setattr(self.options, optkey, os.path.abspath(path))


        # we create a temporal file to store the HMM output
        hmm_output = self.create_temp_file()
        self.press_hmm_library()
        self.run_hmm_scan(hmm_output)
        self.process_hmm_output(hmm_output)
 
    def press_hmm_library(self):
        """Runs ``hmmpress`` on the collection of HMM to make
	it a binary file an make it able to be run with hmm_scan
	"""
        # First, we check whether the binary files exist
        binary_filenames = [self.hmm_file + ".h3" + i for i in ["f", "i", "m", "p"]]
	if not all([os.path.isfile(fi) for fi in binary_filenames]):
            # the binaries do not exist, we create them
            self.log.info("Pressing the HMM library")
            args = []
            args.extend(["-f", " "])
            args.extend([self.hmm_file])

	    cmd = self.get_tool_cmdline("hmmpress", args)
            if not cmd:
                self.log.fatal("cannot find hmmpress in %s" % self.options.hmmpress_path)
                return False
            
            hmmpress_prog = subprocess.Popen(args, stdin=open(os.devnull), stdout=(os.devnull))
            retcode = hmmpress_prog.wait()

            if retcode != 0:
                self.log.fatal("hmmpress exit code was %d, exiting..." % retcode)
                return False

            self.log.info("HMM models pressed")


    def get_tool_cmdline(self, tool_name, args = None):
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
        if not os.path.exists(tool_path):
        # The path does not exist. Assume that this is a file referring to
        # the HMMer tools, but the user has the new one. Try to extract
        # the folder and see if the folder exists.
            base, tool_name = os.path.split(tool_path)

        if os.path.isdir(tool_path):
	# The path exists and points to a folder
	    base = os.path.normpath(tool_path)

        # Okay, first check for the actual tool in the folder
        full_path = os.path.join(tool_path, tool_name)

        if os.path.isfile(full_path):
            return [full_path] + args
        # Nothing succeeded
        return None

    def process_hmm_output(self, hmm_output):
        current_sequence = None
        with open(hmm_output) as f:
	    for line in f:
                if not line.startswith("#"):
                    words = line.split
		    [seq, _, _, model, _, _, score, _, _, _, _, _, _, _, _, sfrom, sto] = words[0:5]

                    if current_sequence is None:
                        current_sequence = seq

    def run_hmm_scan(self, hmm_output_file):
        """Runs ``hmmscan`` on the given sequence file, writing the output to the
        given ``hmm_output_file``.
        Returns ``True`` if the execution was successful, ``False`` otherwise.
        """
        self.log.info("Invoking hmmscan, this might take a long time...")

        args = []
        args.extend(["--cpu", str(self.options.num_threads)])
	args.extend(["--domtblout", hmm_output_file])
        args.extend(["--acc"])
	args.extend(["--incE", str(self.options.evalue)])

        cmd = self.get_tool_cmdline("hmmscan", args)
        if not cmd:
            self.log.fatal("cannot find hmmscan in %s" % self.options.hmmscan_path)
            return False	    

        hmmscan_prog = subprocess.Popen(args, stdin=open(os.devnull), stdout=open(os.devnull))
        retcode = hmmscan_prog.wait()
        if retcode != 0:
            self.log.fatal("hmmscan exit code was %d, exiting..." % retcode)
            return False

        self.log.info("hmmscan returned successfully.")
        return True

if __name__ == "__main__":
    sys.exit(HMMScanApp().run())
