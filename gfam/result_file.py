#!/usr/bin/env python

from gfam.utils import open_anything
import re

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

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

    def write_result_from_dict(self, d, valid_proteins=None, significance=None):
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
                sorted_proteins = sorted_proteins & valid_proteins
            filter_pvalue = False
            if significance is not None:
                filter_pvalue = True
            for protein, goterm_pvalues in valid_proteins.items():
                out.write(protein + "\n")
                for goterm, pvalue in sorted(goterm_pvalues, key=lambda x: x[1], reverse=True):
                    if not filter_pvalue or (filter_pvalue and pvalue < significance):
                        out.write("{}: {}\n".format(goterm, pvalue))
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
                pvalue = pvalue.replace(':', '')
                if pvalue < alpha:
                    list_go_terms.append((goterm, float(pvalue)))
        else:
            if list_go_terms:
                file_content[prev_prot] = list_go_terms
        return file_content

    def get_result_as_dict(self):
        return self._read_result_file(self.file_name)
