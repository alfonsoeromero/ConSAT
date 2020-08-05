import unittest

from gfam.tasks.seqslicer.slice import Slice
from gfam.tasks.seqslicer.slice_builder import SliceBuilder
from gfam.tasks.seqslicer.slice_error import SliceError


class TestSliceBuilder(unittest.TestCase):
    def test_wrong_order_should_raise_exception(self):
        # arrange
        line = "prot_id 2 1"
        
        # act / assert
        with self.assertRaises(SliceError):
            SliceBuilder.from_string(line)

    def test_wrong_number_of_fields_should_raise_exception(self):
        # arrange
        line = ["prot_id 2",
                "prot_id 1 2 3"]
        
        # act / assert
        with self.assertRaises(SliceError):
            SliceBuilder.from_string(line[0])
            SliceBuilder.from_string(line[1])
            
    def test_slice_should_be_perfectly_built(self):
        # arrange
        line = "prot_id 1 2"
        expected_slice = Slice("prot_id", 1, 2)
        
        # act
        obtained_slice = SliceBuilder.from_string(line)
        
        # assert
        self.assertEqual(expected_slice, obtained_slice)
