"""Unit tests for the symbol table."""
import unittest
from unittest.mock import Mock

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import SequenceStore, SliceMapping, SequenceAnnotation
from supergsl.core.types.position import Slice, Position
from supergsl.core.constants import STRAND_CRICK


class SequenceStoreTestCase(unittest.TestCase):
    """Testcases to evaluate the SymbolTable class."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.sequence_store = self.fixtures.sequence_store

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

        self.assertEqual(entry1.sequence_annotations(), [annotation1])


    def test_retrieve_annotations_from_slice(self):
        """"""
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

        results = entry1.sequence_annotations_for_slice(
            Slice.from_five_prime_indexes(40,210))

        expected_results = [annotations[1], annotations[2]]
        self.assertEqual(results, expected_results)
