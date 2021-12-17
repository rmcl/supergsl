"""Tests for the JSON output module."""
from unittest import TestCase
from io import StringIO
from Bio.Seq import Seq

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin.output.genbank_output import GenBankOutput

class GenBankOutputTestCase(TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

        self.compiler_settings = {}
        self.output = GenBankOutput(
            self.fixtures.mk_provider_config(self.compiler_settings))

    def test_build_seq_record_for_assembly(self):
        """Test that the SeqRecord is created correctly."""
        assembly = self.fixtures.mk_assembly(num_parts=2)

        record = self.output.build_seq_record_for_assembly(assembly)

        self.assertEqual(record.seq, assembly.sequence)
        self.assertEqual(record.name, assembly.identifier)
        self.assertEqual(
            record.description,
            'This is a great assembly named asm1')

        self.assertEqual(len(record.features), 2)
        part0_feature = record.features[0]
        self.assertEqual(part0_feature.id, 'part-000')
        self.assertEqual(part0_feature.type, 'part')

        part1_feature = record.features[1]
        self.assertEqual(part1_feature.id, 'part-001')
        self.assertEqual(part1_feature.type, 'part')
