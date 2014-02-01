#!/usr/bin/env python
# -*- coding: utf-8 -*- 
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
import urllib
import tempfile
import shutil
import os
import io
import sys
import math
import re
import operator
import unicodedata
import string
from os import listdir
from collections import Counter, defaultdict
from gfam.scripts import CommandLineApp
from gfam.architecture import ArchitectureFileReaderPerArch as ArchReader

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"


class GetText(CommandLineApp):
    """\
    Usage: %prog [options]

    Application that retrieves the terms associated to a set of sequences
    from PubMed and the description itself of the proteins. It makes a
    weighting operation in order to assign an importance to each one of
    the terms. Given that weight, it would be possible to perform a ranking
    of the terms if needed, or simply to get a vector of textual features for
    each sequence.
    """

    short_name = "get_text"

    my_stopwords = """abundant additional ago analyses analysis analyzed 
                  approximately assigned based
                  carry causing chromosome chromosomes closely common commonly
                  comparative compared comparison complete comprehensive
                  comprehensively confirmed correlated data dataset decay derive
                  derived details detected describe determined differences 
                  divergence divergences diverse diversity draft early encoding 
                  estimated evidence evolution evolutionary existence explain 
                  families family find found fragment full gene genes 
                  genetic genome genomes genomic group groups high highly 
                  identification 
                  identified identifying important infer initial improved
                  includes 
                  including independent independently indicating 
                  interestingly involved isolated isolates laboratory 
                  large largely larger leading level levels lineage 
                  lineages located long loss lost low made
                  maintain major member method methods model molecular modern 
                  number observed obtain occur occurred orfs perspective 
                  perspectives phylogenetic predicted present 
                  previously protein proteins proteome provide putative 
                  recent reflect related remains remaining remarkable 
                  remarkably report represents results revealed reveals sequence 
                  sequenced sequences sequencing set shown similar 
                  single species specific spite strain strains strongly 
                  study studies taxa the type uncharacterized unique validated 
                  wide widely widespread years yields
                    """.split()

    def __init__(self, *args, **kwds):
        super(GetText, self).__init__(*args, **kwds)

    def create_parser(self):
        """Creates the command line parser for this application"""
        parser = super(GetText, self).create_parser()
        parser.add_option("-s", "--stopwords", dest="stopwords_file",
                          metavar="FILE", config_key="file.stopwords",
                          help="Stopwords file. If none is provided," +
                               "a standard one will be used")
        parser.add_option("-S", "--sequences",
                          dest="sequences_file", metavar="FILE",
                          config_key="analysis:find_unassigned/sequences_file",
                          help="FASTA file containing all the "
                               "sequences to analyze")
        parser.add_option("-m", "--mapping", dest="mapping", metavar="FILE",
                          config_key="file.idmapping",
                          help="'idmapping_selected.tab' file, which can " +
                               "be downloaded from " +
                               "ftp://ftp.uniprot.org/pub/databases/uniprot" +
                               "/current_release/knowledgebase/idmapping/")
        parser.add_option("-r", "--seq-id-regexp", metavar="REGEXP",
                          config_key="sequence_id_regexp",
                          help="remap sequence IDs using REGEXP (only " +
                               "used in the FASTA file)")
        parser.add_option("-a", "--arch-file", metavar="FILE", dest="arch_file",
                          help="table with architectures produced by GFam" +
                               " to give final results per architecture")                          
        parser.add_option("-f", "--rdf-file", dest="rdf_file", metavar="FILE",
                          config_key="file.rdffile",
                          help="RDF file containing publication abstract " +
                               "which can be downloaded from http://www." +
                               "uniprot.org/citations/?query=citedin%3a(*) " +
                               "&format=*&compress=yes")
        parser.add_option("-o", "--output-dir", dest="output", metavar="DIR",
                          config_key="folder.output",
                          help="directory where the output files are to be " +
                               "written.")
        parser.add_option("-l", "--lexicon", dest="lexicon_file", metavar="FILE",
                          config_key="file.lexicon", help="file with a previously " +
                          "generated lexicon used to maintain compatibility")
        parser.add_option("-p", "--pubmed-cache", dest="pubmed_cache", metavar="DIR",
                         config_key="folder.pubmed_cache",
                         help="A directory where a cache of PubMed downloaded abstacts " +
                         "will be stored")
        return parser

    def read_stopwords(self, fileName):
        """Reads the stopwords file to a set and returns it.
        If no file is provided, the standard SMART stopword
        list is donwnloaded"""
        if not fileName or not os.path.exists(fileName):
            # we download the SMART list
            url = "http://jmlr.csail.mit.edu/papers/volume5/lewis04a/a11-smart-stop-list/english.stop"
            fd, fileName = tempfile.mkstemp()
            urllib.urlretrieve(url, fileName)

        self.stopwords = set([line.strip() for line in open(fileName)])
        self.stopwords = self.stopwords | set(GetText.my_stopwords)

    def check_not_exists(self, fileName):
        """Checks that a file does not exists. If it does, exists and
        shows an error message"""
        if os.path.exists(fileName):
            self.error("The file " + fileName +
                       " already exists. Cannot overwrite")

    def read_lexicon_file(self, file_name):
        """ Reads a previously generated lexicon file and maintains compatibility
            with its identifiers. The document frequency, however, is discarded
        """
        if not os.path.exists(file_name):
            self.error("The given lexicon file " + file_name + " does not exist")
        for line in open(file_name, "r"):
            word, word_id, _ = line.split()
            self.old_lexicon[word] = int(word_id)

    def run_real(self):
        """Runs the applications"""
        # 1.- we check the compulsory arguments
        if not self.options.sequences_file:
            self.error("the Fasta file should be provided")

        if not self.options.rdf_file:
            self.error("the RDF file should be provided")

        if not self.options.mapping:
            self.error("the Mapping file should be provided")

        self.old_lexicon = dict()
        if self.options.lexicon_file:
            self.read_lexicon_file(self.options.lexicon_file)

        # 2.- we read the stopwords data
        self.read_stopwords(self.options.stopwords_file)

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
        #self.check_not_exists(self.lexicon_file)
        self.freq_file = os.path.join(output_dir, "freq_file_per_sequence")
        #self.check_not_exists(self.freq_file)
        self.weight_file = os.path.join(output_dir, "weight_file_per_sequence")
        #self.check_not_exists(self.weight_file)
        self.ids_file = os.path.join(output_dir, "seq_ids_file")
        #self.check_not_exists(self.ids_file)
        self.text_file = os.path.join(output_dir, "text_file_per_sequence")
        #self.check_not_exists(self.text_file)
        #self.weight_file_arch = os.path.join(output_dir, "weight_file_per_arch")
        #self.check_not_exists(self.weight_file_arch)

        # 4.- we process the rdf file to a cache directory
        if not self.options.pubmed_cache:
            self.cachepubmed = os.path.join(output_dir, "cache_pubmed")
        else:
            self.cachepubmed = self.options.pubmed_cache
        self.log.info("Preprocessing PubMed abstracts")
        if not os.path.exists(self.cachepubmed):
            os.makedirs(self.cachepubmed)
        self.process_rdf_file()

        # 5.- we call the main function to extract text of the protein
        #self.cacheuniprot = os.path.join(output_dir, "cache_uniprot")
        #self.log.info("Extracting protein descriptions")

        #if not os.path.exists(self.cacheuniprot):
        #     os.makedirs(self.cacheuniprot)
        #self.process_fasta_file()

        # 5.1.- the text file is created
        self.log.info("Creating file with all text pieces for each protein")
        self._digits = re.compile('\d')
        self.punctuation = set(string.punctuation) - set(['-'])
        #self._numeral = re.compile(u"^[\+\<\>\~\~\≥\*\=\-\±]*\d+(\.\d*)?([\+\<\>\~\~\≥\*\=\-\±]+(\d+(\.\d*))?)?([nmcμ]m|%|μgml|[mtc]|[kmκ]b|[x×]|umg|mua|fold)*$")
        self._numeral = re.compile(u"^([\+\<\>\~\~\≥\*\=\-\±]*\d+)(\.\d*)?([\+\<\>\~\~\≥\*\=\-\±·]+\d*(\.\d*)?)?([nmcμ]m|%|μgml|[mtc]|[kmκ]b|[x×]|umg|mua|fold)*$")
        self.create_text_file()

        # 6.- once the text has been extracted, the file is preprocessed
        #  indexed and all the output files written
        self.log.info("Indexing the text file")
        self.index_text_file()
        self.log.info("Creating weight file")
        self.weight_freq_file()

        if self.options.arch_file:
            self.log.info("Creating weight arch file")
            self.weight_arch_file()

        # END of the module

    def create_text_file(self):
        seq2pmid = defaultdict(set)
        self.log.info("Reading mapping file")
        pattern = re.compile(";\s*")
        with open(self.options.mapping, "r") as f:
            for line in f:
                fields = line.split("\t")
                if fields[16].strip():
                    seq_id, pmids = fields[0], fields[16]
                    seq2pmid[seq_id] = set(map(int, pattern.split(pmids)))

        self.log.info("Writing text file")
        with io.open(self.text_file, "w", encoding='utf8') as out, open(self.options.sequences_file, 'r') as f:
            cache_pmid = dict()
            for line in f:
                if line.startswith(">"):
                    id_prot = line.split()[0].replace(">", "")
                    if "|" in id_prot:
                        id_prot = id_prot.split("|")[1]
                    text_prot = self.normalize(line.split(" ", 1)[1].split("OS=")[0].strip())
                    
                    text_pubmed = []
                    if id_prot in seq2pmid:
                        for pmid in seq2pmid[id_prot]:
                            if pmid in cache_pmid:
                                text_pubmed.append(cache_pmid[pmid])
                                continue
                            pmid_fi = os.path.join(self.cachepubmed, str(pmid))
                            try:
                                text_pmid = self.normalize(open(pmid_fi, 'r').read())
                                text_pubmed.append(text_pmid)
                                cache_pmid[pmid] = text_pmid
                            except IOError:
                                pass
                    out.write(u"{0} {1} {2}\n".format(id_prot, text_prot, " ".join(text_pubmed)))

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
                if "<title>" in line:
                    opentitle = "</title>" not in line
                    if not opentitle:
                        title = line.split("<title>")[1].split("</title>")[0]
                    else:
                        title = line.split("<title>")[1]
                elif "http://www.ncbi.nlm.nih.gov/pubmed/" in line:
                    pubmedid = line.split("pubmed/")[1].split("\"")[0]
                elif "</rdf:Description>" in line:
                    # we process here the record
                    filename = os.path.join(self.cachepubmed, pubmedid)
                    with open(filename, "w") as out:
                        out.write("{} {}".format(title, " ".join(abstract)))
                    pubmedid, title, abstract = "", "", []
                    opentitle, openabstract = False, False
                elif "<rdfs:comment>" in line:
                    if "</rdfs:comment>" in line:
                        abstract.append(line.split("<rdfs:comment>")[1].
                                        split("</rdfs:comment>")[0])
                    else:
                        abstract.append(line.split("<rdfs:comment>")[1])
                        openabstract = True
                elif opentitle:
                    if "</title>" not in line:
                        title += line
                    else:
                        title += line.split("</title>")[0]
                elif openabstract:
                    if "</rdfs:comment>" not in line:
                        abstract.append(line)
                    else:
                        abstract.append(line.split("</rdfs:comment>")[0])
                        openabstract = False

    def process_fasta_file(self):
        with open(self.options.sequences_file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    id_prot = line.split()[0].replace(">", "")
                    txt_protein = line.split(" ", 1)[1].split("OS=")[0].strip()
                    filename = os.path.join(self.cacheuniprot, id_prot)
                    with open(filename, "w") as out:
                        out.write("{}".format(txt_protein))

    def normalize(self, s):
        """This is a pre-processing which is made to any text
        before it is written to the disk (and prior to be
        indexed). Unicode data is normalized, the stopwords
        (and small tokens) are removed as well as numbers and URLs, 
        and it is set to lowercase.
        """
        u_s = ''.join((c for c in unicodedata.normalize('NFD', unicode(s.lower(), "UTF-8"))
            if unicodedata.category(c) != 'Mn'))
        u_s = ''.join([c if not c in self.punctuation else ' ' for c in u_s])
        return ' '.join([x.strip('-') for x in u_s.split() if len(x) > 2 and len(x) < 30 and 
            x not in self.stopwords
            and not self.is_number(x) 
            and not self.is_numeral(x)
            and not x.startswith('http') 
        ])

    def tokenize(self, s):
        """Tokenizes and preprocesses a certain string into
        a list of strings. Yes, it can (should be) improved"""
        u_s = unicode(s, "UTF-8")
        return u_s.split()

    def is_numeral(self, s):
        """Checks if the string is sort of a number (i.e. an interval, or a percentage)
        """
        if not self._digits.search(s):
            return False

        return bool(self._numeral.match(s))

    def is_number(self, s):
        """Checks if a string is a number"""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def index_text_file(self):
        freq = defaultdict(int)
        if len(self.old_lexicon) == 0:
            word_id = dict()
        else:
            word_id = self.old_lexicon
        seq_ids = []
        # open freq file as output
        with open(self.freq_file, 'w') as output:
            for line in open(self.text_file, 'r'):
                [seq_id, rest] = line.split(" ", 1)
                seq_ids.append(seq_id)

                line = [seq_id]
                for token, frequency in Counter(self.tokenize(rest)).iteritems():
                    tok_id = 0
                    if token not in word_id:
                        tok_id = len(word_id)
                        word_id[token] = tok_id
                    else:
                        tok_id = word_id[token]
                    line.append("{}:{}".format(tok_id, frequency))
                    freq[tok_id] += 1
                output.write("{}\n".format(" ".join(line)))

        self.N = len(seq_ids)

        # write lexicon
        with io.open(self.lexicon_file, 'w', encoding='utf8') as lexicon:
            for word, id_word in word_id.iteritems():
                if id_word in freq:
                    lexicon.write(u"{0} {1} {2}\n".format(word, id_word, freq[id_word]))

        # write seq ids file
        with open(self.ids_file, 'w') as proteins_file:
            for protein in seq_ids:
                proteins_file.write("{}\n".format(protein))

    def weight_arch_file(self):
        # 1.- read architectures to a dictionary of sets (key=architecture,
        # value set of sequences
        self.log.info("Inverting dict")
        seqs_per_arch = dict((arch, set(prots)) for (arch, prots) in ArchReader(self.options.arch_file))
        arch_per_seq = dict()
        for arch, prots in seqs_per_arch.iteritems():
            arch_per_seq.update([(prot, arch) for prot in prots])

        # 2.- main loop
        self.log.info("Creating weight arch file")
        vec_per_arch = dict((arch, defaultdict(float)) for arch in seqs_per_arch.keys())
        for line in open(self.weight_file, 'r'):
            fields = line.split()
            if not fields[0] in arch_per_seq or len(fields) == 0:
                continue
            arch = arch_per_seq[fields[0]]
            vec = vec_per_arch[arch] 
                
            for field in fields[1:]:
                term, weight = field.split(':')
                vec[int(term)] += float(weight)

        self.log.info("Printing weight arch file")
        for arch in seqs_per_arch: 
            print arch, " ", " ".join(["{0}:{1:.5f}".format(t,w) for t,w in vec_per_arch[arch].iteritems()])

    def weight_freq_file(self):
        """Reads the frequency file and produces the frequency file"""
        # we read the lexicon file
        freq = dict()
        for line in open(self.lexicon_file, 'r'):
            [_, numid, doc_freq] = line.split()
            freq[int(numid)] = int(doc_freq)

        # we read the frequency file and write it
        with open(self.weight_file, 'w') as output_file:
            idf = dict([(numid, math.log(float(self.N)/float(df))) for numid, df in freq.iteritems()])
            for line in open(self.freq_file, 'r'):
                if " " in line.strip():
                    seq_id, rest = line.split(" ", 1)
                    features = [(int(f[0]), float(f[1])) for f in [i.split(":") for i in rest.split()]]
                    output_file.write("{} {}\n".format(seq_id, " ".join(["{}:{:.5f}".format(feat[0], feat[1]*idf[feat[0]]) for feat in features])))
                else:
                    output_file.write("{}\n".format(seq_id))
        self.log.info("Weight file done!")

    def check_not_exists(self, filename):
        try:
            with open(filename):
                self.error("The file " + filename + " already exists")
        except IOError:
            pass

if __name__ == "__main__":
    sys.exit(GetText().run())
