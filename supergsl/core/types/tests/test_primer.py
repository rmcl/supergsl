"""Tests for the output core module."""
import unittest
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.types.primer import PrimerPair


class BuiltinTypeTestCase(unittest.TestCase):
    """Test the behavior of SuperGSL built in types."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.symbol_table = self.fixtures.mk_global_symbol_table()
        self.sequence_store = self.symbol_table.lookup('sequences')

    def test_primer_pair_from_sequence(self):
        """Test that PrimerPair type can be initialized from a pair of sequences."""

        forward_entry = self.sequence_store.add_from_reference(Seq('ATGCCCCCCCATAGA'))
        reverse_entry = self.sequence_store.add_from_reference(Seq('TTAGACACATGGGAC'))

        result = PrimerPair.from_sequence_entries(forward_entry, reverse_entry)
        self.assertEqual(result.forward.sequence, forward_entry.sequence)
        self.assertEqual(result.reverse.sequence, reverse_entry.sequence)
