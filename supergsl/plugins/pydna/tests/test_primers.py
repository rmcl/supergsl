"""Tests for the output core module."""
import unittest
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.pydna.primers import ExtractionPrimerBuilder

class ExtractionPrimerBuilderTestCase(unittest.TestCase):
    """Test the behavior of PyDNA based primer builder"""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.primer_builder = ExtractionPrimerBuilder(
            sequence_store=self.fixtures.sequence_store)

    def test_build_primers_for_part(self):
        """Test that primers are created for a part."""

        _, part = self.fixtures.mk_part('AwesomePart', 200, False)

        self.assertEqual(part.extraction_primers, None)

        amplicon = self.primer_builder.build_primers_for_part(part)

        primer_pair = part.get_extraction_primers()
        self.assertEqual(
            primer_pair.forward.sequence,
            amplicon.forward_primer.seq)
        self.assertEqual(
            primer_pair.reverse.sequence,
            amplicon.reverse_primer.seq)
