import re

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"

import itertools

class IndexedReadOnlyFile(object):
    """ We provide a mechanism to easily access a set of blocks from
        a text file. We define a block as a set of lines starting with
        a string, matching a certain `regex`, and ending either with the
        following block, or with the EOF.

        Once initialized, the method `get_keys()` is available, returning
        the set of valid keys. To access the list of lines from an object
        "obj" of this class, we can now make obj[key]
    """
    def __init__(self, file_name, regex):
        self.file_name = file_name
        self.regex = re.compile(regex)
        self.key_to_offset = self.__preload_cache()

    def get_keys(self):
        """ Return the set of keys which are accessible from this file.
        """
        return self.key_to_offset.keys()

    def __getitem__(self, k):
        """ Returns the list of items (file lines) associated with
            a certain key. The line corresponding to the key is not
            returned (you should call `get_line_for_key` instead).
            If the key is not a valid key an empty line is returned.
        """
        if k in self.key_to_offset:
            with open(self.file_name, 'r') as input:
                input.seek(self.key_to_offset[k], 0)
                output = []
                for line in input:
                    l = line.strip()
                    if l:
                        output.append(l)
                    else:
                        break
                #l = itertools.takewhile(lambda line : line.strip(), iter(input))
                return output
        else: 
            return []

    def get_line_for_key(self, key):
        """ Returns the whole line corresponding to a certain key. 
            This may be useful if the key is a subset of the line
            and we want to retrieve extra information from the it.
        """
        if key in self.key_to_offset:
            return line
        return ""

    def __preload_cache(self):
        with open(self.file_name, 'r') as input:
            return dict((line.rstrip(), input.tell()) 
                    for line in iter(input.readline, '') 
                        if self.regex.match(line))
#
#
