#!/usr/bin/env python
"""All-against-all BLASTing of a given set of sequences.

"""

from __future__ import with_statement

import os
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
    Usage: %prog [options] sequences_file

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
	hmm_output = create_temp_file()
	press_hmm_library()
        run_hmm_scan(hmm_output)
	process_hmm_output(hmm_output)
 

    def press_hmm_library(self):
        pass

    def process_hmm_output(self, hmm_output):
        pass

    def run_hmm_scan(self, hmm_output_file):
        """Runs ``hmmscan`` on the given sequence file, writing the output to the
	given ``hmm_output_file``.

        Returns ``True`` if the execution was successful, ``False`` otherwise.
        """
        self.log.info("Invoking hmmscan, this might take a long time...")

        # where the unused output of hmmscan will be printed
	# this file will be ignored
	dummy_output = create_temp_file()

        args = []
        args.extend(["--cpu", str(self.options.num_threads)])
	args.extend(["-o", dummy_output])
	args.extend(["--tblout", hmm_output_file])
	args.extend(["--incE", str(self.options.evalue)])
	args.extend(["--noali"])

#####################################################################################################

        args = self.get_blast_cmdline("", args)
        if not args:
            self.log.fatal("cannot find blastall in %s" % self.options.blastall_path)
            return False

#####################################################################################################

        hmmscan_prog = subprocess.Popen(args, stdin=open(os.devnull))
        retcode = hmmscan_prog.wait()
        if retcode != 0:
            self.log.fatal("hmmscan exit code was %d, exiting..." % retcode)
            return False

        self.log.info("hmmscan returned successfully.")
        return True

if __name__ == "__main__":
    sys.exit(HMMScanApp().run())
