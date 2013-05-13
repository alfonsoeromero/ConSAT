#!/usr/bin/env python
"""Application that retrieves a set of weighted terms for the sequences
in order to have clues of what is their function. The sources where the
text is downloaded from are PubMed abstracts and the text appearing in 
the description of the protein itself.

After being retrieved, the text is written to an intermediate file which
is then parsed and preprocessed (case folding, punctuation symbols removed, 
stopwords removal, etc). A simple tf-idf is then performed for each set
of terms.

The terms associated to a certain sequence can be easily retrieved with the
files produced.
"""

from collections import defaultdict

import urllib
import tempfile, shutil, os

import sys
import math
from os import listdir

from gfam.scripts import CommandLineApp

from collections import Counter

__author__  = "Alfonso E. Romero"
__email__   = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

class GetText(CommandLineApp):
    """\
    Usage: %prog [options]

    Application that retrieves the terms associated to a set of sequences
    from PubMed and the description itself of the proteins. It makes a weighting
    operation in order to assign an importance to each one of the terms. Given that
    weight, it would be possible to perform a ranking of the terms if needed, or
    simply to get a vector of textual features for each sequence.
    """

    short_name = "get_text"

    def __init__(self, *args, **kwds):
        super(GetText).__init__(*args, **kwds)

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(GetText, self).create_parser()
        parser.add_option("-s", "--stopwords", dest="stopwords_file",
                metavar="FILE",help="Stopwords file. If none is provided,"+
                                    "a standard one will be used")
        parser.add_option("-S", "--sequences",
                dest="sequences_file", metavar="FILE",
                help="FASTA file containing all the sequences to analyze")
        parser.add_option("-m", "--mapping", dest="mapping", metavar="FILE",
                help="'idmapping_selected.tab' file, which can be downloaded " +
                "from ftp://ftp.uniprot.org/pub/databases/uniprot/" +
                "current_release/knowledgebase/idmapping/")
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                help="remap sequence IDs using REGEXP (only used in the "+
                     "FASTA file)")
        parser.add_option("-d", "--rdf-file", dest="rdf_file", metavar="FILE",
                help="RDF file containing publication abstract which can be "+
                "downloaded from http://www.uniprot.org/citations/"+
                "?query=citedin%3a(*)&format=*&compress=yes")
        parser.add_option("-o", "--output-dir", dest="output", metavar="DIR",
                help="directory where the output files are to be written")
        return parser


    def read_stopwords(self, fileName):
        """Reads the stopwords file to a set and returns it.
        If no file is provided, the standard SMART stopword
        list is donwnloaded"""
        if not os.path.exists(fileName):
            # we download the SMART list
            url = "http://jmlr.csail.mit.edu/papers/volume5/lewis04a/a11-smart-stop-list/english.stop"
            fd, fileName = tempfile.mkstemp()
            urllib.request.urlretrieve(url, fileName)

        self.stopwords = set([line.strip() for line in open(fileName)])


    def check_not_exists(self, fileName):
        """Checks that a file does not exists. If it does, exists and
        shows an error message"""
        if os.path.exists(fileName):
            self.error("The file " + fileName + " already exists. Cannot overwrite")


    def run_real(self):
        """Runs the applications"""
        
        # 1.- we check the compulsory arguments
        if not self.options.sequences_file :
            self.error("the Fasta file should be provided")

        if not self.options.rdf_file:
            self.error("the RDF file should be provided")

        if not self.options.mapping:
            self.error("the Mapping file should be provided")

        # 2.- we read the stopwords data
        self.read_stopwords(self.stopwords_file)

        # 3.- we check that the output directory is writable, and create the
        # names of the output files
        if not self.options.output:
            self.error("No output directory has been provided")
        if not os.path.exists(self.options.output):
            os.makedirs(self.options.output)
        elif not os.path.isdir(self.options.output):
            self.error("The route " + self.options.output + 
                       " exists, but it is not a directory")

        # we setup the output file names
        output_dir = os.path.normpath(self.options.output)
        self.lexicon_file = os.path.join(output_dir, "lexicon")
        check_not_exists(self.lexicon_file)
        self.freq_file = os.path.join(output_dir, "freq_file")
        check_not_exists(self.freq_file)
        self.weight_file = os.path.join(output_dir, "weight_file")
        check_not_exists(self.weight_file)
        self.ids_file = os.path.join(output_dir, "seq_ids_file")
        check_not_exists(self.ids_file)
        self.text_file = os.path.join(output_dir, "text_file")
        check_not_exists(self.text_file)

        # 4.- we process the rdf file to a cache directory
        self.cachepubmed = os.path.join(output_dir, "cache_pubmed")
        self.process_rdf_file()

        # 5.- we call the main function to extract text of the protein
        self.cacheuniprot = os.path.join(output_dir, "cache_uniprot")
        self.process_fasta_file()

        # 5.1.- the text file is created
        self.create_text_file()

        # 6.- once the text has been extracted, the file is preprocessed
        #  indexed and all the output files written
        self.index_text_file()
        self.weight_freq_file()
        # END of the module

    def create_text_file(self):
        seq2pmid = defaultdict(set)
        with open(self.options.mapping, "r") as f:
            for line in f:
                fields = line.split("\t")
                if len(fields[16].strip()) > 0:
                    seq_id, pmids = fields[0], fields[16]
                    seq2pmid[seq_id] = set(pmids.split(";\s*"))

        with open(self.text_file, "w") as out:
            for id_prot in listdir(self.cacheuniprot):
                file_prot = os.path.join(self.cacheuniprot, id_prot)
                text_prot = open(file_prot, 'r').read()
                text_pubmed = []
                if id_prot in seq2pmid:
                    for pmid in seq2pmid[id_prot]:
                        pmid_file = op.path.join(self.cachepubmed, pmid)
                        try:
                            text_pmid = open(pmid_file, 'r').read()
                            text_pubmed.append(text_pmid)
                        except IOError:
                            pass
                out.write("%s " % id_prot)
                out.write("%s " % text_prot)
                out.write("%s\n" % " ".join(text_pubmed))
        
    def process_rdf_file(self):
        """Process the rdf file, extracting the relevant fields into
        individual text files.
        This parser assumes that the XML file is relatively well formed
        (this might change in the future and some error tolerance has
        been added)
        """
        with open(self.options.rdf_file, "r") as f:
            opentitle = False
            openabstract = False
            pubmedid, title, abstract = "", "", []

            for line in f:
                if line.find("<title>") != -1:
                    opentitle = line.find("</title>") == -1
                    if not opentitle:
                        title = line.split("<title>")[1].split("</title>")[0]
                    else:
                        title = line.split("<title>")[1]
                elif line.find("http://www.ncbi.nlm.nih.gov/pubmed/") != -1:
                    pubmedid = line.split("pubmed/")[1].split("\"")[0]
                elif line.find("</rdf:Description>") != -1:
                    # we process here the record
                    filename = os.path.join(self.cachepubmed, pubmedid)
                    with open(filename, "w") as out:
                        out.write("%s " % title)
                        out.write("%s" % " ".join(abstract))
                    pubmedid, title, abstract = "", "", []
                    opentitle, openabstract = False, False
                elif line.find("<rdfs:comment>") != -1:
                    if line.find("</rdfs:comment>") != -1:
                        abstract.append(line.split("<rdfs:comment>")[1].split("</rdfs:comment>")[0])
                    else:
                        abstract.append(line.split("<rdfs:comment>")[1])
                        openabstract = True
                elif opentitle:
                    if line.find("</title>") == -1:
                        title += line
                    else:
                        title += line.split("</title>")[0]
                elif openabstract:
                    if line.find("</rdfs:comment>") == -1:
                        abstract.append(line)
                    else:
                        abstract.append(line.split("</rdfs:comment>")[0])
                        openabstract = False

    def process_fasta_file(self):
        with open(self.options.sequences_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    id_prot = line.split()[0].replace(">", "")
                    text_protein = line.split(" ",1)[1].split("OS=")[0].strip()
                    filename = os.path.join(self.cacheuniprot, id_prot)
                    with open(filename, "w") as out:
                        out.write("%s" % text_protein)

    def tokenize(self, string):
        """Tokenizes and preprocesses a certain string into
        a list of strings. Yes, it can (should be) improved"""
        return [x for x in string.translate(None, (",./;'?&()")).lower().split() 
                   if x not in self.stopwords and not is_number(x) and len(x) > 2]

    def is_number(self,s):
        """Checks if a string is a number"""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def index_text_file(self):
        freq = defaultdict(int)
        word_id = dict()
        seq_ids = []

        # open text file
        with open(self.freq_file, 'w') as output:
            for line in open(self.text_file,'r'):
                [seq_id, rest] = line.split(" ", 1)
                seq_ids.append(seq_id)

                num_tokens = []
                for token in tokenize(rest):
                    if token not in word_id:
                        word_id[token] = len(word_id)
                    
                    num_token.append(word_id[token])
                sorted_x = sorted(Counter(num_tokens).iteritems(), key=operator.itemgetter(1))
                for word_id, count in sorted_x:
                    freq[word_id] += 1
                    # write freq file
                    output.write(str(word_id) + ":" + str(count) + " ")
                output.write('\n')                        

        # write lexicon
        with open(self.lexicon_file, 'w') as lexicon:
            for word, id_word in word_id.items():
                lexicon.write(" ".join([word, str(id_word), str(freq[id_word]), "\n"]))
                
        # write seq ids file
        with open(self.ids_file, 'w') as proteins_file:
            for protein in seq_ids:
                lexicon.write(protein + "\n")

    def weight_freq_file(self):
        """Reads the frequency file and produces the frequency file"""
        # we read the lexicon file
        freq = dict()
        for line in open(self.lexicon_file,'r'):
            [numid, _, doc_freq] = text.split()
            freq[int(numid)] = int(doc_freq)

        # we read the frequency file
        with open(self.freq_file, 'w') as output_file:
            for line in open(self.freq_file,'r'):
                fields = line.split()
                seq_id = fields.pop(0)
                features = [i.split(":") for i in fields]
                output_features = [seq_id]
                for feat in features:
                    l = [feat[0], str( math.log( float(self.N)/freq(int(feat[0])) ) * float(feat[1]) )]
                    output_features.append(":".join(l))
                output_file.write(" ".join(output_features))
                output_file.write('\n')

if __name__ == "__main__":
    sys.exit(GetText().run())
    
