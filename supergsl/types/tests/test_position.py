import unittest
import mock
from Bio.Seq import Seq
from supergsl.core.constants import THREE_PRIME, FIVE_PRIME
from supergsl.types.position import SeqPosition

class SeqPositionTestCase(unittest.TestCase):
    """Test case for SeqPosition."""
    maxDiff = None

    def setUp(self):
        self.seq_ex1 = Seq(self.example_seq)


    def test_seq_position_string_representation(self):
        sp = SeqPosition.from_reference(
            x=100,
            rel_to=THREE_PRIME,
            approximate=True,
            reference=self.seq_ex1
        )

        self.assertEquals(str(sp), "3'+100bp Approx:True (Abs Pos: 100)")

    def test_get_absolute_position_in_reference_no_parent_part_three_prime(self):
        sp = SeqPosition.from_reference(
            x=200,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=self.seq_ex1
        )

        ref, pos = sp.get_absolute_position_in_reference()

        self.assertEquals(ref, self.seq_ex1)
        self.assertEquals(pos, 200)

    def test_get_absolute_position_in_reference_no_parent_part_five_prime(self):
        sp = SeqPosition.from_reference(
            x=255,
            rel_to=FIVE_PRIME,
            approximate=False,
            reference=self.seq_ex1
        )

        ref, pos = sp.get_absolute_position_in_reference()

        self.assertEquals(ref, self.seq_ex1)
        self.assertEquals(pos, len(self.seq_ex1) - 255)

    def test_get_absolute_position_in_reference_with_parent_part(self):
        parent_pos = SeqPosition.from_reference(
            x=150,
            rel_to=FIVE_PRIME,
            approximate=False,
            reference=self.seq_ex1
        )

        child_pos = parent_pos.get_relative_position(50)

        ref, pos = child_pos.get_absolute_position_in_reference()
        self.assertEquals(ref, self.seq_ex1)
        self.assertEquals(pos, len(self.seq_ex1) - 100)


    example_seq  = '''AAATAGCCCTCATGTACGTCTCCTCCAAGCCCTGTTGTCTCTTACCCGGATGTTCAACCAAAAGCTACTTACTACCTTTATTTTATGTTTACTTTTTATAGGTTGTCTTTTTATCCCACTTCTTCGCACTTGTCTCTCGCTACTGCCGTGCAACAAACACTAAATCAAAACAATGAAATACTACTACATCAAAACGCATTTTCCCTAGAAAAAAAATTTTCTTACAATATACTATACTACACAATACATAATCACTGACTTTCGTAACAACAATTTCCTTCACTCTCCAACTTCTCTGCTCGAATCTCTACATAGTAATATTATATCAAATCTACCGTCTGGAACATCATCGCTATCCAGCTCTTTGTGAACCGCTACCATCAGCATGTACAGTGGTACCCTCGTGTTATCTGCAGCGAGAACTTCAACGTTTGCCAAATCAAGCCAATGTGGTAACAACCACATCTCCGAAATCTGCTCCAAAAGATATTCCAGTTTCTGCCGAAATGTTTTATTGTAGAACAGCCCTATCAGCATCGACAGGAATGCCGTCCAATGCGGCACTTTAGATGGGGTAACTCCCAGCGCAAGCTGATCTCGCAAGTGCATTCCTAGACTTAATTCATATCTGCTCCTCAACTGTCGATGATGCCTGCTAAACTGCAGCTTGACGTACTGCGGACCCTGCAGTCCAGCGCTCGTCATGGAACGCAAACGCTGAAAAACTCCAACTTTCTCGAGCGCTTCCACAAAGACCGTATCGTCTTTTGCCTCCCATTCTTCCCGGCACTTTTTTTCGTCCCAGTTCAAAAAGTACTGCAGCACCTCTGTCTTCGATTCACGCAAGTTGCTCCATACTTTATAATACAACTCTTTGATCTGCCTTCCAGACATGCGGAAAACTTGGCTCCCTTGCTTGCCTCTTGTCGAATCCAATACACTAATTGTTTCTCTTCTTCTAGTAATGGCCAGGTACCAAGCATAATTTCTCTGTATCTGAGAGTAGATCTCTCTCCTTTTTACGCTAAAATATTTCAAATATCCTACAGGGTCCCCATGATATGGCTCGATGTCTTCCAAGTATTCTTTGTATTCCTCATCATTTCGCAGCATTCTCTCCACAGCTAGTGCTTCCCAAGCTATCCTCCGATACGATACTTTCTGGCCAGCCCAACAGACACAGAGCCCGAACATCTTTTGACAGCCCTTGCATAATCCGTATTGTGTGAATACTCCCTCTGGGCAGAAGTATATGTCAATACCATAGAGGAAAAGATGTTTAATTTCGTCAGACCGAAATCCAAGAAACTGTAAGACATTCATATTCTCGGAAGTATTGGGAAATTGTGCTTTCAGTTTCTTTCTCTCTAGGAAAACCATTTGACTCCCTTTCCGCTTATACGACTCTTTGTTAATGTCGGTGACTGGATGGAATCTATTATCCTCAGCATTGCCATCTTTATTGGCGTCCTCCTTGGCACTAGCGTTGGTACTTTCAGTGGTAGTGGCATTAGTGCTGGAGTTGGTGCTAGCAGTGGTAGTGGCATTAGTGCTGGAGTTGGTGCTAGCAGTGGTAGTAGCACTAGTGTTGGAGTCGGTACTTTCGGTGGTAGTAGCACTAGTGTTGGAGTTGGTACTTTCAGT'''
