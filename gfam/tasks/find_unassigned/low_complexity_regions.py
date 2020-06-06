import re
from collections import defaultdict

from gfam.utils import open_anything


def read_low_complexity_regions(file_name: str,
                                sequence_id_regexp: str = "") ->\
        defaultdict:
    low_complexity_regions = defaultdict(list)
    current_prot_id = ""

    if sequence_id_regexp:
        regexp = re.compile(sequence_id_regexp)
    else:
        regexp = None

    for line in open_anything(file_name):
        line = line.strip()
        if not line:
            continue
        if line[0] == ">":
            current_prot_id = line.split()[0][1:]
            if regexp is not None:
                current_prot_id = regexp.sub(r'\g<id>', current_prot_id)
        else:
            (left, _, right) = line.split()
            left = int(left)
            right = int(right)
            low_complexity_regions[current_prot_id].append((left, right))
    return low_complexity_regions
