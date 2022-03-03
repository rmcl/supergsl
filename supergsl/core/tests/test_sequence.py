"""Unit tests for the symbol table."""
import unittest
from unittest.mock import Mock

from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import (
    SequenceStore,
    SequenceEntry,
    SliceMapping,
    SequenceAnnotation
)
from supergsl.core.types.position import Slice, Position, AbsolutePosition
from supergsl.core.constants import STRAND_CRICK, THREE_PRIME, FIVE_PRIME
from supergsl.core.exception import SequenceStoreError


class SequenceStoreTestCase(unittest.TestCase):
    """Testcases to evaluate the SymbolTable class."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.sequence_store = self.fixtures.sequence_store

    def test_initialize_sequence_entry(self):
        """SequenceEntry must be initialized with reference or parents, but not both."""

        with self.assertRaises(SequenceStoreError) as store_error:
            SequenceEntry(self.fixtures.sequence_store, 123)

        self.assertEqual(
            str(store_error.exception),
            'Must specify either parents or reference.')

        with self.assertRaises(SequenceStoreError) as store_error_2:
            SequenceEntry(
                self.fixtures.sequence_store,
                123,
                reference=Seq('ATGC'),
                parent_links=['hi'])

        self.assertEqual(
            str(store_error_2.exception),
            'Must only specify parents or a reference sequence, but not both')

    def test_add_from_reference(self):
        seq1 = 'ATGAAACACAAATTTAGACACAGAGTAGACATACGATGGAA'
        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1)

        self.assertEqual(entry1, store.lookup(entry1.id))
        self.assertEqual(entry1.sequence, seq1)


    def test_english_text_example(self):
        """Use english sentences to demonstrate the sequence stores concatenate."""
        ref1_str = 'RUSSELL IS A GREAT DUDE YAY FOR HIM LETS TEST CHOPPING UP THIS SEQUENCE'
        ref2_str = 'TODAY IS SATURDAY AND I AM SITTING HERE READING A NOVEL ABOUT HORSES'
        ref3_str = 'JOE CAT IS SITTING ON THE COUCH PURING ABOUT PUREEE AND ENJOYING AN EXPRESSO'
        expected_output = 'RUSSELL IS READING A NOVEL ABOUT HORSES AND ENJOYING AN EXPRESSO'

        store = self.fixtures.sequence_store
        ref1 = store.add_from_reference(ref1_str)
        ref2 = store.add_from_reference(ref2_str)
        ref3 = store.add_from_reference(ref3_str)

        desired_slices = [
            [ref1, [0, 11]],
            [ref2, [40, 68]],
            [ref3, [51, 76]]
        ]

        part_map = []
        trg_pos = 0
        new_trg_pos = 0
        for slice_detail in desired_slices:
            reference_entry = slice_detail[0]
            start_pos = slice_detail[1][0]
            end_pos = slice_detail[1][1]

            new_trg_pos += (end_pos - start_pos)
            src_slice = Slice.from_five_prime_indexes(start_pos, end_pos)
            trg_slice = Slice.from_five_prime_indexes(trg_pos, new_trg_pos)

            trg_pos = new_trg_pos

            part_map.append(SliceMapping(
                reference_entry,
                src_slice,
                trg_slice
            ))

        new_entry = store.concatenate(part_map)

        self.assertEqual(new_entry.sequence, expected_output)

    def test_slice_part_with_reference(self):
        """Simple slice behavior of a parent sequence with only a reference sequence."""

        reference_sequence = self.fixtures.mk_random_dna_sequence(500)
        sequence_entry = self.fixtures.mk_sequence_entry(reference_sequence)

        new_entry = self.sequence_store.slice(
            sequence_entry,
            Slice.from_five_prime_indexes(100,150)
        )

        self.assertEqual(new_entry.sequence, reference_sequence[100:150])

    def test_slice_grandparent_part(self):
        """Slice behavior on a part that has itself been sliced from a part with a reference sequence.

        P1:     |---X-------------X-------------|
        P2:         |---X----X----|
        P3:             |----|

        """
        reference_sequence = self.fixtures.mk_random_dna_sequence(500)
        grandparent_entry = self.fixtures.mk_sequence_entry(reference_sequence)

        parent_entry = self.sequence_store.slice(
            grandparent_entry,
            Slice.from_five_prime_indexes(100,250)
        )

        child_entry = self.sequence_store.slice(
            parent_entry,
            Slice.from_five_prime_indexes(50,100)
        )

        self.assertEqual(parent_entry.sequence, reference_sequence[100:250])
        self.assertEqual(child_entry.sequence, reference_sequence[150:200])

    def test_slice_reverse_strand(self):
        """Test slice behavior when desired slice is on CRICK strand."""

        reference_sequence = self.fixtures.mk_random_dna_sequence(500)
        parent_entry = self.fixtures.mk_sequence_entry(reference_sequence)

        child_entry = self.sequence_store.slice(
            parent_entry,
            Slice.from_five_prime_indexes(100,250, strand=STRAND_CRICK)
        )

        self.assertEqual(child_entry.sequence, reference_sequence.reverse_complement()[100:250])

    def test_slice_grandparent_part_with_reverse_strand(self):
        """Slice a part that has itself been sliced and resides on reverse strand of reference sequence.

        P1:     |---X-------------X-------------|
        P1RC:   |-------------X-------------X---| <-- reverse complement of parent
        P2:                   |---X----X----|
        P3:                       |----|

        """
        reference_sequence = self.fixtures.mk_random_dna_sequence(500)
        grandparent_entry = self.fixtures.mk_sequence_entry(reference_sequence)

        parent_entry = self.sequence_store.slice(
            grandparent_entry,
            Slice.from_five_prime_indexes(100,250, strand=STRAND_CRICK)
        )

        child_entry = self.sequence_store.slice(
            parent_entry,
            Slice.from_five_prime_indexes(50,100)
        )

        expected_parent_sequence = reference_sequence.reverse_complement()[100:250]
        self.assertEqual(parent_entry.sequence, expected_parent_sequence)

        expected_child_sequence = expected_parent_sequence[50:100]
        self.assertEqual(child_entry.sequence, expected_child_sequence)


class SequenceAnnotationTestCase(unittest.TestCase):
    """Testcases to evaluate the handling of sequence annotations."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.sequence_store = self.fixtures.sequence_store

    def test_add_and_retrieve_annotation(self):
        """Create an annotation, add it to a reference sequence and then retrieve it."""
        seq1 = 'ATGAAACACAAATTTAGACACAGAGTAGACATACGATGGAA'

        annotation1 = SequenceAnnotation.from_five_prime_indexes(0,20, ['HELLO'], {
            'payload': 'party'
        })

        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1, annotations=[annotation1])

        self.assertEqual(entry1.annotations(), [annotation1])

    def test_retrieve_annotations_from_slice(self):
        """Return only annotations that are found within a slice."""
        seq1 = self.fixtures.mk_random_dna_sequence(2000)
        annotations = [
            SequenceAnnotation.from_five_prime_indexes(0,20, ['HELLO'], {}),
            SequenceAnnotation.from_five_prime_indexes(50,200, ['YO'], {}),
            SequenceAnnotation.from_five_prime_indexes(70,190, ['YO2'], {}),
            SequenceAnnotation.from_five_prime_indexes(55,215, ['YO3'], {}),
            SequenceAnnotation.from_five_prime_indexes(500,1000, ['YO4'], {})
        ]

        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1, annotations=annotations)

        results = entry1.annotations_for_slice(
            Slice.from_five_prime_indexes(40,210))

        expected_results = [annotations[1], annotations[2]]
        self.assertEqual(results, expected_results)

    def test_annotations_from_parent_parts(self):
        """Return annotations that have been defined on a parent part."""
        seq1 = self.fixtures.mk_random_dna_sequence(2000)
        annotations_on_parent = [
            SequenceAnnotation.from_five_prime_indexes(0, 20, ['HELLO'], {}),
            SequenceAnnotation.from_five_prime_indexes(50, 200, ['YO'], {}),
            SequenceAnnotation.from_five_prime_indexes(70, 190, ['YO2'], {}),
        ]
        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1, annotations=annotations_on_parent)

        part_slice = Slice.from_five_prime_indexes(50, 250, strand=STRAND_CRICK)
        annotations_on_child = [
            SequenceAnnotation.from_five_prime_indexes(10, 30, ['X1'], {}),
            SequenceAnnotation.from_five_prime_indexes(80, 190, ['X2'], {}),
        ]
        entry2 = store.slice(entry1, part_slice, annotations=annotations_on_child)

        parent_start = AbsolutePosition(entry2.sequence_length, 50, False)
        expected_results = [
            annotations_on_parent[1].derive_from_absolute_start_position(parent_start),
            annotations_on_child[0],
            annotations_on_parent[2].derive_from_absolute_start_position(parent_start),
            annotations_on_child[1]
        ]
        results = entry2.annotations()
        print(results)
        self.assertEqual(results, expected_results)

    def test_annotations_from_parent_part_reverse_strand(self):
        """Test that we return annotations from parent parts on reverse strand."""

        seq1 = self.fixtures.mk_random_dna_sequence(2000)
        a1_slice = Slice(
            Position(-20, THREE_PRIME, False),
            Position(0, THREE_PRIME, False),
            strand=STRAND_CRICK)
        annotations_on_parent = [
            SequenceAnnotation(a1_slice, ['HELLO-FROM-REV-STRAND'], {}),
        ]
        store = self.fixtures.sequence_store
        entry1 = store.add_from_reference(seq1, annotations=annotations_on_parent)
        entry2 = store.slice(
            entry1,
            Slice.from_five_prime_indexes(0, 250))

        parent_start = AbsolutePosition(entry2.sequence_length, 0, False)
        results = entry2.annotations()
        self.assertEqual(results, [
            annotations_on_parent[0].derive_from_absolute_start_position(parent_start)
        ])

    def test_entry_link_sequence(self):
        """A entry link sequence should return the correct sequence from the parent entry."""

        reference_sequence = self.fixtures.mk_random_dna_sequence(500)
        sequence_entry = self.fixtures.mk_sequence_entry(reference_sequence)

        new_entry = self.sequence_store.slice(
            sequence_entry,
            Slice.from_five_prime_indexes(100,150)
        )

        self.assertEqual(
            new_entry.parent_links[0].sequence,
            reference_sequence[100:150])
