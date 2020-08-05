import unittest

from gfam.sequence import Sequence, SeqRecord


class TestSequence(unittest.TestCase):
    def test_sequence_should_be_a_string(self):
        # arrange
        s = "ACCTAAC"

        # act
        sut = Sequence(s)

        # assert
        self.assertEqual(sut, s)


class TestSeqRecord(unittest.TestCase):
    def test_seqrecord_should_be_creatable_with_default_values(self):
        # arrange
        se = "ACCTAAC"
        description = "description"
        id = "id"
        name = "name"

        # arrange
        s = SeqRecord(se)
        s1 = SeqRecord(se, id)
        s2 = SeqRecord(se, id, name)
        s3 = SeqRecord(se, id, name, description)

        # assert
        self.assertEqual(s.seq, se)
        self.assertEqual(s.seq, s1.seq)
        self.assertEqual(s1.seq, s2.seq)
        self.assertEqual(s3.seq, s2.seq)

        self.assertEqual(s1.id, id)
        self.assertEqual(s1.id, s2.id)
        self.assertEqual(s2.id, s3.id)

        self.assertEqual(s2.name, name)
        self.assertEqual(s2.name, s3.name)

        self.assertEqual(s3.description, description)
        self.assertEqual(s3.description, description)

    def test_fragment_should_produce_a_fragment_with_the_entire_sequence(self):
        # arrange
        se = "ACCTAAC"
        id = "my_prot"
        len_prot = len(se)

        # act
        s = SeqRecord(se, id)
        frag = s.fragment()

        # assert
        self.assertEqual(s.seq, frag.seq)
        self.assertEqual(frag.id, f"{id}:1-{len_prot}")
