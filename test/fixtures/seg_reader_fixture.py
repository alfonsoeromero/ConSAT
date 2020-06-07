import os
import re
from collections import defaultdict
from test.fixtures.find_unassigned_fixtures import FindUnassignedFixture
from typing import Dict, Set, Tuple


class SEGReaderFixture:
    """Class used to manually read SEG regions
        into a dict[protein, List[Tuple[int, int]]]
    """

    def __init__(self):
        current_dir = os.path.dirname(__file__)
        fixture = FindUnassignedFixture()
        self.data_dir = os.path.join(current_dir, os.pardir, "data")
        self.seg_file = fixture.get_low_complexity_regions_file()

    def get_dict_low_complexity_regions(self) ->\
            Dict[str, Set[Tuple[int, int]]]:
        low_complexity_regions = defaultdict(set)
        current_prot_id = ""

        sequence_id_regexp = r"(\w+\|)(?P<id>\w+)(\|\w+)+"
        regexp = re.compile(sequence_id_regexp)

        for line in open(self.seg_file):
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                current_prot_id = line.split()[0][1:]
                current_prot_id = regexp.sub(r'\g<id>', current_prot_id)
            else:
                (left, _, right) = line.split()
                left = int(left)
                right = int(right)
                low_complexity_regions[current_prot_id].add((left, right))
        return low_complexity_regions
