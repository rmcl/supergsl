"""Test SuperGSL builtin types."""
from typing import Union, Optional
from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures

from supergsl.core.types import SuperGSLType
from supergsl.core.types.position import Slice
from supergsl.core.types.builtin import (
    Collection,
    SliceAndInvertCollection,
    SliceInvertMixin,
    NucleotideSequence,
    AminoAcidSequence
)


class SliceablePartMock(SuperGSLType, SliceInvertMixin):
    """Mock class for a Sliceable type."""
    def __init__(self):
        self.slice_called = False
        self.last_slice = None
        self.invert_called = False

    def slice(
        self,
        part_slice : Union[Slice, str],
        identifier : Optional[str] = None
    ) -> 'SuperGSLType':
        """Return a new part representing a sliced region."""
        self.slice_called = True
        self.last_slice = part_slice
        return self

    def invert(self, identifier : Optional[str] = None):
        """Return a new Part with sequence on reverse strand."""
        self.invert_called = True
        return self


class BuiltinBasicTypesTestCase(TestCase):
    """Tests for AminoAcidSequence and NucleotideSequence"""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()


    def test_amino_acid_sequence_repr(self):
        """Test the repr for a AminoAcidSequence."""
        aa_seq = self.fixtures.mk_random_amino_acid_sequence(525)
        aa_seq_entry = self.fixtures.mk_sequence_entry(aa_seq)
        amino_acid_seq = AminoAcidSequence(aa_seq_entry)

        self.assertEqual(
            repr(amino_acid_seq),
            f'AminoAcidSequence: {aa_seq_entry.sequence[0:55]}...{aa_seq_entry.sequence[-3:]}')

    def test_nucleotide_sequence_repr(self):
        """Test the repr from a NucleotideSequence."""
        seq_entry = self.fixtures.mk_random_dna_sequence_entry(500)
        nucleotide_sequence = NucleotideSequence(seq_entry)

        self.assertEqual(
            repr(nucleotide_sequence),
            f'NucleotideSequence: {seq_entry.sequence[0:55]}...{seq_entry.sequence[-3:]}')




class CollectionTestCase(TestCase):
    """TestCase for the Collection type."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_collection_len(self):
        """len() returns the correct number of items in collection."""
        collection1 = self.fixtures.mk_part_collection(num_parts=7)
        self.assertEqual(len(collection1), 7)

    def test_collection_repr(self):
        """Test that the output of the collection repr is expected."""
        collection1 = self.fixtures.mk_part_collection(num_parts=7)

        expected_output = f'Collection (count: {len(collection1._items)})\n'
        for item_idx, item in enumerate(collection1):
            expected_output += f'  {item_idx}. {item} \n'

        self.assertEqual(str(collection1), expected_output)

    def test_iterate_over_collection(self):
        """Test collection interator."""
        collection1 = self.fixtures.mk_part_collection()
        self.assertEqual(list(collection1), collection1._items)


class SliceAndInvertCollectionTestCase(TestCase):
    """TestCase for the SliceAndInvertCollection type."""

    def test_iterate_attempts_to_slice_and_invert(self):
        """Validate that the iter method slices and inverts each part in the collection."""

        collection1 = Collection([
            SliceablePartMock(),
            SliceablePartMock()
        ])
        s1 = Slice.from_entire_sequence()

        slice_inv_collection = SliceAndInvertCollection(collection1, s1, True)

        for item in slice_inv_collection:
            self.assertEqual(item.slice_called, True)
            self.assertEqual(item.last_slice, s1)
            self.assertEqual(item.invert_called, True)
