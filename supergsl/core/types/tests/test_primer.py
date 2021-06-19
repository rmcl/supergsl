"""Tests for the output core module."""
import unittest
from Bio.Seq import Seq
from supergsl.core.types.primer import PrimerPair


class BuiltinTypeTestCase(unittest.TestCase):
    """Test the behavior of SuperGSL built in types."""
    maxDiff = None

    def test_primer_pair_from_sequence(self):
        """Test that PrimerPair type can be initialized from a pair of sequences."""

        forward = Seq('ATGCCCCCCCATAGA')
        reverse = Seq('TTAGACACATGGGAC')

        result = PrimerPair.from_sequences(forward, reverse)
        self.assertEqual(result.forward.sequence, forward)
        self.assertEqual(result.reverse.sequence, reverse)
