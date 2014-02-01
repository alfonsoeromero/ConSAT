#!/usr/bin/env python

from gfam.utils import open_anything
from gfam.indexed_file import IndexedReadOnlyFile
import re
import math 

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

try:
    from scipy.stats import chisqprob
except ImportError:
    def chisqprob(chi, df):
        """Return prob(chisq >= chi, with df degrees of
            freedom).
            df must be even.
            copypasted from 
            http://www.linuxjournal.com/files/linuxjournal.com/linuxjournal/articles/064/6467/6467s2.html
        """
        assert df & 1 == 0
        # XXX If chi is very large, exp(-m) will underflow to 0.
        m = chi / 2.0
        sum = term = math.exp(-m)
        for i in range(1, df//2):
            term *= m / i
            sum += term
        # With small chi and large df, accumulated
        # roundoff error, plus error in
        # the platform exp(), can cause this to spill
        # a few ULP above 1.0. For
        # example, chi2P(100, 300) on my box
        # has sum == 1.0 + 2.0**-52 at this
        # point.  Returning a value even a teensy
        # bit over 1.0 is no good.
        return min(sum, 1.0)

class ResultFileFilter(object):
    """Filters a result file by p-value, keeping only 
        those entries under a certain significance value
    """
    def __init__(self, input_file):
        self.input_file = input_file

    def filter(self, out_file, confidence=0.05):
        reader = ResultFileReader(self.input_file, confidence)
        writer = ResultFileWriter(out_file)
        for key in reader.get_keys():
            l = reader[key]
            names = reader.get_go_names()
            writer.write_single_register(key, l, significance=confidence, go_names=names)
        writer.close()

class ResultFileCombiner(object):
    """Combines two result files using Fisher's method
    """
    def __init__(self, file_name1, file_name2):
        self.file_name1 = file_name1
        self.file_name2 = file_name2

    def combine(self, out_file, confidence=0.05):
        infile1 = ResultFileReader(self.file_name1)
        infile2 = ResultFileReader(self.file_name2)
        keys1 = set(infile1.get_keys())
        keys2 = set(infile2.get_keys())
        all_keys = keys1 | keys2
        out = ResultFileWriter(out_file)
        names = dict()
        for key in sorted(all_keys):
            l_out = []
            if key in keys1 and key in keys2:
                p1 = dict(infile1[key])
                p2 = dict(infile2[key])
                for goterm in set(p1.keys()) | set(p2.keys()):
                    if goterm in p1 and goterm in p2:
                        l_out.append((goterm, self.combine_fisher(p1[goterm], p2[goterm])))
                    elif goterm in p1:
                        l_out.append((goterm, p1[goterm]))
                    else:
                        l_out.append((goterm, p2[goterm]))
            elif key in keys1:
                l_out = infile1[key]
            else:
                l_out = infile2[key]
            names.update(infile1.get_go_names())
            names.update(infile2.get_go_names())
            out.write_single_register(key, l_out, significance=confidence, go_names=names)
        out.close()

    def combine_fisher(self, pvalue1, pvalue2):
        if pvalue1 == 0.0 or pvalue2 == 0.0:
            return 0.0
        else:
            chi = -2.0 * (math.log(pvalue1) + math.log(pvalue2))
            p_out = chisqprob(chi, 4)
            return p_out

class ResultFileWriter(object):
    """Class implementing a result file writer, with the
    format:
    protein_id
        p-value: GO:XXXXXXX
        p-value: GO:XXXXXXX
        ...

    protein_id2
        p-value: ...
        ...
    
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.out = open(self.file_name, "w")

    def close(self):
        """ Closes the output stream where the results have been
            written
        """
        self.out.close()

    def write_single_register(self, key, l, significance=None, go_names=None):
        """ Writes a single register represented by a `key` and a list `l` of
            pairs (goterm, pvalue). We can add a `significance` cutoff and
            a dictionary with the names of the GO terms
        """
        filter_pvalue = False or significance is not None
        self.out.write("{}\n".format(key))
        for goterm, pvalue in sorted(l, key=lambda x: x[1], reverse=True):
            if not filter_pvalue or (filter_pvalue and pvalue < significance):
                if go_names is not None and goterm in go_names:
                    self.out.write("\t{:.5f}: {} ({})\n".format(pvalue, goterm, go_names[goterm]))
                else:
                    self.out.write("\t{:.5f}: {}\n".format(pvalue, goterm))
        self.out.write("\n")            

    def write_result_from_dict(self, d, valid_proteins=None, significance=None, go_names=None):
        """Writes a result file from a dictionary structure where,
        for each protein identifier we get a list of tuples 
        (goterm, pvalue) representing the assignment itself.
        There is no need for up-propagation as the result files
        should be already up-propagated. If we specify a `significance`
        value, then only entries with a p-value less than `significance`
        will be considered.
        If `valid_proteins` is specified, only sequences in the 
        intersection of that list and the keys of the dictionary will
        be written. If we specify a `significance`
        value, then only entries with a p-value less than `significance`
        will be considered.
        """
        sorted_proteins = sorted(d.keys())
        if valid_proteins is not None:
            sorted_proteins = sorted(set(sorted_proteins) & set(valid_proteins))
        filter_pvalue = False
        if significance is not None:
            filter_pvalue = True
        for protein in sorted_proteins:
            goterm_pvalues = d[protein]
            self.write_single_register(protein, goterm_pvalues, significance, go_names)

class ResultFileReader(object):
    """Class implementing a parser of the
    GFam result file, i.e., a text file with the
    following format:

    protein_id
        p-value: GO:XXXXXXX (function)
        p-value: GO:XXXXXXX (function)
        ...

    protein_id2
        p-value: ...
        ...
    """
    def __init__(self, file_name, significance=None):
        """ Prepares a ResultFileReader to read such from a `file_name`.
            If we specify a certain `significance`, we will only consider
            those contents with a pvalue less than this `significance`.
        """
        self.file = IndexedReadOnlyFile(file_name, "^\w+")
        self.keys = self.file.get_keys()
        self.pval_regex = re.compile("^\d")
        self.go_names = dict()
        if significance is not None:
            self.alpha = significance
        else:
            self.alpha = 1000.0

    def get_keys(self):
        """ Return the set of protein-ids of the file
        """
        return self.keys

    def get_result_as_dict(self):
        """ Retrieves the whole dataset as a dictionary.
            Not recommended if the file is too large.
        """
        d = dict()
        for key in self.keys:
            d[key] = self.__getitem__(key)
        return d

    def __getitem__(self, key):
        """ Gets the set of GO terms and p-values for a
            certain key (which should be a real key in
            the file).
        """
        list_go_terms = []
        for line in self.file[key]:
            l = line.strip()
            if self.pval_regex.match(l):
                pvalue, goterm = l.split(' ')[0:2]
                pvalue = float(pvalue.replace(':', ''))
                if pvalue < self.alpha:
                    list_go_terms.append((goterm, pvalue))
                    if goterm not in self.go_names:
                        name = l.split('(', 1)[1][0:-1]
                        self.go_names[goterm] = name
        return list_go_terms

    def get_go_names(self):
        if not self.go_names:
            for key in self.keys:
                self.__getitem__(key)            
        return self.go_names

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print "ERROR: arguments should be p-value_file1 p-value_file2 output_combined output_filtered"
        sys.exit(-1)

    in1, in2, out1, out2 = sys.argv[1:]
    print "Reading first file"
    reader1 = ResultFileReader(in1, 0.04)
    for key in reader1.get_keys():
        reader1[key]
    print "read ", len(reader1.get_go_names()), " goterms"

    print "Reading second file"
    reader2 = ResultFileReader(in2)
    reader2.get_result_as_dict()
    print "read ", len(reader2.get_go_names()), " goterms"

    print "Combining files and writing it"
    comb = ResultFileCombiner(in1, in2)
    comb.combine(out1)

    print "Filtering combined file"
    filterer = ResultFileFilter(out1)
    filterer.filter(out2, confidence=0.05)

    print "Everything done!"
