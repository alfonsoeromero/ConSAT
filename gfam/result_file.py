#!/usr/bin/env python

""" Utilities and functions to process a result
    file (i.e. triplets protein identifier, GO term
    and p-values)
"""

import math
import re

from gfam.indexed_file import IndexedReadOnlyFile

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

try:
    from scipy.stats import chisqprob
except ImportError:
    def chisqprob(chi, deg):
        """Return prob(chisq >= chi, with deg degrees of
            freedom).
            deg must be even.
            copypasted from
            http://www.linuxjournal.com/files/linuxjournal.com/
                   linuxjournal/articles/064/6467/6467s2.html
        """
        assert deg & 1 == 0
        # XXX If chi is very large, exp(-m) will underflow to 0.
        half_chi = chi / 2.0
        summation = term = math.exp(-half_chi)
        for i in range(1, deg//2):
            term *= half_chi / i
            summation += term
        # With small chi and large deg, accumulated
        # roundoff error, plus error in
        # the platform exp(), can cause this to spill
        # a few ULP above 1.0. For
        # example, chi2P(100, 300) on my box
        # has summation == 1.0 + 2.0**-52 at this
        # point.  Returning a value even a teensy
        # bit over 1.0 is no good.
        return min(summation, 1.0)


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
            li = reader[key]
            names = reader.get_go_names()
            writer.write_single_register(key, li,
                                         significance=confidence,
                                         go_names=names)
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
                pvalues_1 = dict(infile1[key])
                pvalues_2 = dict(infile2[key])
                for goterm in set(pvalues_1.keys()) | set(pvalues_2.keys()):
                    if goterm in pvalues_1 and goterm in pvalues_2:
                        l_out.append((goterm,
                                      self.combine_fisher(pvalues_1[goterm],
                                                          pvalues_2[goterm])))
                    elif goterm in pvalues_1:
                        l_out.append((goterm, pvalues_1[goterm]))
                    else:
                        l_out.append((goterm, pvalues_2[goterm]))
            elif key in keys1:
                l_out = infile1[key]
            else:
                l_out = infile2[key]
            names.update(infile1.get_go_names())
            names.update(infile2.get_go_names())
            out.write_single_register(key, l_out, significance=confidence,
                                      go_names=names)
        out.close()

    def combine_fisher(self, pvalue1, pvalue2):
        """ Combine two p-values using Fihser's method. See
            https://en.wikipedia.org/wiki/Fisher%27s_method for
            more details
        """
        if pvalue1 == 0.0 or pvalue2 == 0.0:
            return 0.0
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

    def write_single_register(self, key, l_r,
                              significance=None, go_names=None):
        """ Writes a single register represented by a `key` and a list `l_r` of
            pairs (goterm, pvalue). We can add a `significance` cutoff and
            a dictionary with the names of the GO terms
        """
        self.out.write("{}\n".format(key))
        if significance is not None:
            sorted_l = sorted([y for y in l_r if y[1] < significance],
                              key=lambda x: x[1])
        else:
            sorted_l = sorted(l_r, key=lambda x: x[1])

        for goterm, pvalue in sorted_l:
            if go_names is not None and goterm in go_names:
                self.out.write("\t{:.5f}: {} ({})\n".format(pvalue,
                                                            goterm,
                                                            go_names[goterm]))
            else:
                self.out.write("\t{:.5f}: {}\n".format(pvalue, goterm))
        self.out.write("\n")


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
        self.file = IndexedReadOnlyFile(file_name, r"^\w+")
        self.keys = self.file.get_keys()
        self.pval_regex = re.compile(r"^\d")
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
        results = dict()
        for key in self.keys:
            results[key] = self.__getitem__(key)
        return results

    def __getitem__(self, key):
        """ Gets the set of GO terms and p-values for a
            certain key (which should be a real key in
            the file).
        """
        list_go_terms = []
        for line in [li for li in self.file[key] if self.pval_regex.match(li)]:
            pvalue, goterm = line.split(' ', 2)[0:2]
            pvalue = float(pvalue[0:-1])
            if pvalue < self.alpha:
                list_go_terms.append((goterm, pvalue))
                if goterm not in self.go_names:
                    name = line.split('(', 1)[1][0:-1]
                    self.go_names[goterm] = name
        return list_go_terms

    def get_go_names(self):
        if not self.go_names:
            for key in self.keys:
                self.__getitem__(key)
        return self.go_names
