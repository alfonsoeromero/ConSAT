#!/usr/bin/env python

from gfam.utils import open_anything
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
        reader = ResultFileReader(self.input_file)
        data = reader.get_result_as_dict(confidence)
        names = reader.get_go_names()
        writer = ResultFileWriter(out_file)
        writer.write_result_from_dict(data, significance=confidence, go_names=names)

class ResultFileCombiner(object):
    """Combines two result files using Fisher's method
    """
    def __init__(self, file_name1, file_name2):
        self.file_name1 = file_name1
        self.file_name2 = file_name2

    def combine(self, out_file, confidence=0.05):
        infile1 = ResultFileReader(self.file_name1)
        infile2 = ResultFileReader(self.file_name2)
        d1 = infile1.get_result_as_dict()
        d2 = infile2.get_result_as_dict()
        all_proteins = set(d1.keys()) | set(d2.keys())
        names1 = infile1.get_go_names()
        names2 = infile2.get_go_names()
        names = dict(names1.items() + names2.items())
        dout = dict()
        for prot in sorted(all_proteins):
            if prot in d1 and prot in d2:
                l_out = []
                p1 = dict(d1[prot])
                p2 = dict(d2[prot])
                for goterm in set(p1.keys()) | set(p2.keys()):
                    if goterm in p1 and goterm in p2:
                        l_out.append((goterm, self.combine_fisher(p1[goterm], p2[goterm])))
                    elif goterm in p1:
                        l_out.append((goterm, p1[goterm]))
                    else:
                        l_out.append((goterm, p2[goterm]))
                dout[prot] = l_out                
            elif prot in d1:
                dout[prot] = d1[prot]
            else:
                dout[prot] = d2[prot]
        out = ResultFileWriter(out_file)
        out.write_result_from_dict(dout, significance=confidence, go_names=names)

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
        with open(self.file_name, "w") as out:
            sorted_proteins = sorted(d.keys())
            if valid_proteins is not None:
                sorted_proteins = sorted(set(sorted_proteins) & set(valid_proteins))
            filter_pvalue = False
            if significance is not None:
                filter_pvalue = True
            for protein in sorted_proteins:
                goterm_pvalues = d[protein]
                out.write(protein + "\n")
                for goterm, pvalue in sorted(goterm_pvalues, key=lambda x: x[1], reverse=True):
                    if not filter_pvalue or (filter_pvalue and pvalue < significance):
                        if go_names is not None and goterm in go_names:
                            out.write("\t{:.5f}: {} ({})\n".format(pvalue, goterm, go_names[goterm]))
                        else:
                            out.write("\t{:.5f}: {}\n".format(pvalue, goterm))
                out.write("\n")            

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
    def __init__(self, file_name):
        self.file_name = file_name

    def _read_result_file(self, file_name, significance=None):
        """Reads the result file into a dictionary structure where,
        for each protein identifier we get a list of tuples 
        (goterm, pvalue) representing the assignment itself.
        There is no need for up-propagation as the result files
        should be already up-propagated. If we specify a `significance`
        value, then only entries with a p-value less than `significance`
        will be considered.
        """
        file_content = dict()
        self.go_names = dict()
        list_go_terms = []
        prev_prot = ""
        prot_regex = re.compile("^[A-Z]") 
        pval_regex = re.compile("^\d")
        if significance is not None:
            alpha = significance
        else:
            alpha = 1000.0

        for line in open_anything(file_name):
            l = line.strip()
            if prot_regex.match(l):
                if len(prev_prot) > 0:
                    file_content[prev_prot] = list_go_terms
                prev_prot = l
                list_go_terms = []
            elif pval_regex.match(l):
                pvalue, goterm = l.split(' ')[0:2]
                pvalue = float(pvalue.replace(':', ''))
                if pvalue < alpha:
                    list_go_terms.append((goterm, pvalue))
                    if goterm not in self.go_names:
                        name = l.split('(', 1)[1][0:-1]
                        self.go_names[goterm] = name
        else:
            if list_go_terms:
                file_content[prev_prot] = list_go_terms
        return file_content

    def get_result_as_dict(self, significance=None):
        return self._read_result_file(self.file_name, significance)

    def get_go_names(self):
        return self.go_names

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print "ERROR: arguments should be p-value_file1 p-value_file2 output_combined output_filtered"
        sys.exit(-1)

    in1, in2, out1, out2 = sys.argv[1:]
    print "Reading first file"
    reader1 = ResultFileReader(in1)
    reader1.get_result_as_dict()
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
