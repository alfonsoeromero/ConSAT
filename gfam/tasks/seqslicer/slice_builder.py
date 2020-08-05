from gfam.tasks.seqslicer.slice import Slice
from gfam.tasks.seqslicer.slice_error import SliceError


class SliceBuilder:
    @staticmethod
    def from_string(line: str) -> Slice:
        fields = line.strip().split()
        seq_id = fields[0]
        left: int = 1
        right: int = -1
        num_fields = len(fields)

        if num_fields <= 2 or num_fields > 3:
            raise SliceError(
                f"Slice found: '{line}' with invalid number of fields")
        else:
            left, right = map(int, fields[1:])
            if left >= right:
                raise SliceError(
                    f"Slice found: '{line}' should have 'right' > 'left'")

        return Slice(seq_id, left, right)
