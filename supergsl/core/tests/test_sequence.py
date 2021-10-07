"""Unit tests for the symbol table."""
import unittest
from unittest.mock import Mock

from supergsl.core.sequence import SequenceStore
from supergsl.core.types.slice import Slice, Position


class SequenceStoreTestCase(unittest.TestCase):
    """Testcases to evaluate the SymbolTable class."""

    def test_add_from_reference(self):
        store = SequenceStore()
        seq1 = 'ATGAAACACAAATTTAGACACAGAGTAGACATACGATGGAA'
        entry1 = store.add_from_reference(seq1)

        self.assertEqual(entry1, store.lookup(entry1.id))
        self.assertEqual(entry1.sequence, seq1)


    def test_english_text_example(self):
        """Use english sentences to demonstrate the sequence stores concatenate."""
        ref1_str = 'RUSSELL IS A GREAT DUDE YAY FOR HIM LETS TEST CHOPPING UP THIS SEQUENCE'
        ref2_str = 'TODAY IS SATURDAY AND I AM SITTING HERE READING A NOVEL ABOUT HORSES'
        ref3_str = 'JOE CAT IS SITTING ON THE COUCH PURING ABOUT PUREEE AND ENJOYING AN EXPRESSO'
        expected_output = 'RUSSELL IS READING A NOVEL ABOUT HORSES AND ENJOYING AN EXPRESSO'

        store = SequenceStore()
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

            src_slice = Slice(Position(start_pos), Position(end_pos))

            new_trg_pos += (end_pos - start_pos)
            trg_slice = Slice(Position(trg_pos), Position(new_trg_pos))
            trg_pos = new_trg_pos

            part_map.append((
                reference_entry,
                src_slice,
                trg_slice
            ))

        new_entry = store.concatenate(part_map)

        self.assertEqual(new_entry.sequence, expected_output)
