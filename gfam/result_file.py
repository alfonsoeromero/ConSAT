#!/usr/bin/env python

from gfam.utils import open_anything
import re

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

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

    def _read_result_file(self, file_name):
        """Reads the result file into a dictionary structure where,
        for each protein identifier we get a list of tuples 
        (goterm, pvalue) representing the assignment itself.
        There is no need for up-propagation as the result files
        should be already up-propagated
        """
        file_content = dict()
        list_go_terms = []
        prev_prot = ""
        prot_regex = re.compile("^[A-Z]") 
        pval_regex = re.compile("^\d")

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
                list_go_terms.append((goterm, float(pvalue)))
        else:
            if list_go_terms:
                file_content[prev_prot] = list_go_terms
        return file_content

    def get_result_as_dict(self):
        return self._read_result_file(self.file_name)
