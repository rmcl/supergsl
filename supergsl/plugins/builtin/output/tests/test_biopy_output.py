"""Tests for the JSON output module."""
from unittest import TestCase
from io import StringIO
from Bio.Seq import Seq

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import SequenceAnnotation
from supergsl.core.types.role import PROMOTER, TERMINATOR
from supergsl.plugins.builtin.output.biopy import (
    SeqRecordAssemblyOutput,
    create_seq_record_from_sequence_entry,
    build_seq_record_for_assembly
)

class BioPythonSeqRecordOutputTestCase(TestCase):
    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_create_seq_record_from_sequence_entry(self):
        """Create a SeqRecord from the given SequenceEntry."""

        seq1 = self.fixtures.mk_random_dna_sequence(2000)

        annotations_on_parent = [
            SequenceAnnotation.from_five_prime_indexes(0, 20, [PROMOTER], {
                'label': 'promoter',
            }),
            SequenceAnnotation.from_five_prime_indexes(70, 190, [TERMINATOR], {}),
        ]
        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1, annotations=annotations_on_parent)

        seq_record = create_seq_record_from_sequence_entry(entry1, 'AWESOME')

        self.assertEqual(seq_record.name, 'AWESOME')
        self.assertEqual(str(seq_record.seq), str(entry1.sequence))

        promoter_feature = seq_record.features[0]
        self.assertEqual(promoter_feature.strand, 1)
        self.assertEqual(promoter_feature.location.start, 0)
        self.assertEqual(promoter_feature.location.end, 20)
        self.assertEqual(promoter_feature.type, 'promoter')

        terminator_feature = seq_record.features[1]
        self.assertEqual(terminator_feature.strand, 1)
        self.assertEqual(terminator_feature.location.start, 70)
        self.assertEqual(terminator_feature.location.end, 190)
        self.assertEqual(terminator_feature.type, 'terminator')


class SeqRecordAssemblyOutputTestCase(TestCase):
    """Test the behavior of GenBank Assembly output"""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

        self.compiler_settings = {}
        self.output = SeqRecordAssemblyOutput(
            self.fixtures.mk_provider_config(self.compiler_settings))

    def test_build_seq_record_for_assembly(self):
        """Test that the SeqRecord is created correctly."""
        assembly = self.fixtures.mk_assembly(num_parts=2)

        record = build_seq_record_for_assembly(assembly)

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
