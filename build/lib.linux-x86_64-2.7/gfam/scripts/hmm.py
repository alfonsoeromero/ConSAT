#!/usr/bin/env python
"""HMM application"""

import optparse
import sys
import os
import errno
from subprocess import call

from gfam import fasta
from collections import defaultdict
from gfam.scripts import CommandLineApp
from gfam.utils import open_anything
from gfam.sequence import SeqRecord


__author__  = "Alfonso E. Romero"
__email__   = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2012, Alfonso E. Romero"
__license__ = "GPL"

class HMM(CommandLineApp):
    """\
    Usage: %prog [options] [sequences_file] [cluster_table_file]

    Learn a set of Hidden Markov Models over each set of sequences
    representing a cluster found in previous stages.

    sequences_file is the file with all the sequences involved in the 
    clustering, in Fasta format

    cluster_table_file is the set of clusters, one by line in a tab-delimited
    file, where the first field is the cluster name, and the rest of fields
    are the sequences names

    """

    short_name = "hmm"

    def __init__(self, *args, **kwds):
        super(HMM, self).__init__(*args, **kwds)

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(HMM, self).create_parser()
        parser.add_option("--dir",
                dest="directory", metavar="FILE",
                help="directory where the HMMs will be stored",
                config_key="generated/file.dir_hmms",
                default=None)

        return parser

    def read_table(self, table_file):
        table = defaultdict()

        table_file_desc = open(table_file, 'r')
        for line in table_file_desc:
            fields = line.strip().split("\t")
            table[fields[0]] = fields[1:]
        table_file_desc.close()

        return table 

    def create_directories(self):
        # 1.- we create the output directory
        maindir = self.options.directory
        try:
            os.makedirs(maindir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        if not maindir.endswith(os.path.sep):
            maindir += os.path.sep

        # 2.- inside it, we create the sequences directory
        self.sequences_dir = maindir + "sequences"
        try:
            os.makedirs(self.sequences_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        # 3.- ... the alignments directory
        self.alignments_dir = maindir + "alignments"
        try:
            os.makedirs(self.alignments_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        # 4.- ... and the hmms directory
        self.hmms_dir = maindir + "hmms"
        try:
            os.makedirs(self.hmms_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def separate_sequences(self, table, sequences_file):
        """Separates the fasta sequence file in individual fasta files,
           one per cluster"""
        reader = fasta.Parser(open_anything(sequences_file))
        seqs = dict(((seq.id, seq) for seq in reader))

        for cluster_name,cluster_seqs in table.iteritems():
            output_file_name = self.sequences_dir + os.path.sep + cluster_name
            output_fd = open(output_file_name + ".faa", "w")
            writer = fasta.Writer(output_fd)

            for sequence in cluster_seqs:
                obj = SeqRecord(seqs[sequence].seq, sequence, "", "")
                writer.write(obj)

            output_fd.close()

    def building_alignments(self, cluster_names):
        """Builds the alignments, one for each cluster name"""
        seq_dir = self.sequences_dir + os.path.sep
        align_dir = self.alignments_dir + os.path.sep

        for cluster in cluster_names:
            seq = seq_dir + cluster + ".faa"
            align = align_dir + cluster + ".sto"
            #print >> sys.stderr, "CLUSTER: ", cluster
            call(["clustalo", "--infile=" + seq, "--outfmt=st", "--outfile=" + align])
	# no errors are checked, in the future the call should be checked

    def building_hmms(self, cluster_names):
        """Builds the HMMs, one for each alignment"""             
        hmm_dir = self.hmms_dir                  
        align_dir = self.alignments_dir
        self.hmm_files = []
            
        for cluster in cluster_names:
            hmm = os.path.join(hmm_dir, cluster) + ".hmm"
            align = os.path.join(align_dir, cluster) + ".sto"

            call(["hmmbuild", hmm, align], \
                    stdout=open(os.devnull, 'wb')) 
            self.hmm_files.append(hmm)
        # TODO: no errors are checked, in the future the call should be checked

    def join_hmms(self):
        """Joins alltogether the produces HMM files into a single one.
        As the HMM files are text-based, the union is basically a
        concatenation of the individual files
        """
        all_models_file = os.path.join(self.options.directory, "all.hmm")
        with open(all_models_file, "w") as out:
            for hmm_file in self.hmm_files:
                for line in open(hmm_file, "r"):
                    out.write(line)

    def run_real(self):
        """Runs the HMM application"""
        if len(self.args) != 2:
            self.error("expected exactly two input file names")

        sequences_file, table_file = self.args

        # TODO: require the presence of hmmbuild and clustalo

        self.log.info("Reading cluster table from %s..." % table_file)
        table = self.read_table(table_file)

        self.log.info("Creating directories...")
        self.create_directories()

        self.log.info("Separating FASTA in one file per cluster...")
        self.separate_sequences(table, sequences_file)

        self.log.info("Building alignments...")
        self.building_alignments(table.keys())

        self.log.info("Building HMMs")
        self.building_hmms(table.keys())

        self.log.info("Joining HMM files into single one")
        self.join_hmms()

if __name__ == "__main__":
    sys.exit(HMM().run())
