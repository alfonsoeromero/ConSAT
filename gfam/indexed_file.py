import re

__author__ = "Alfonso E. Romero"
__email__ = "aeromero@cs.rhul.ac.uk"
__copyright__ = "Copyright (c) 2013, Alfonso E. Romero"
__license__ = "GPL"


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
            with open(self.file_name, 'r') as input_f:
                input_f.seek(self.key_to_offset[k], 0)
                output = []
                for line in input_f:
                    line_strip = line.strip()
                    if line_strip:
                        output.append(line_strip)
                    else:
                        break
                return output
        else:
            return []

    def __preload_cache(self):
        """ Preloads the internal cache
        """
        with open(self.file_name, 'r') as input_f:
            return dict((line.rstrip(), input_f.tell())
                        for line in iter(input_f.readline, '')
                        if self.regex.match(line))
