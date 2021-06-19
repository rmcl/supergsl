import unittest
from Bio.Seq import Seq
from supergsl.core.ast import (
    Slice as AstSlice,
    SlicePosition as AstSlicePosition,
    SymbolReference as AstSymbolReference
)
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.constants import (
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)

from supergsl.core.parts.slice import convert_slice_position_to_seq_position


class ResolvePartSlicePassTestCase(unittest.TestCase):
    """Test case for `ResolvePartSlicePass`."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_convert_slice_position_to_seq_position_with_relative_to_start(self):
        reference, part = self.fixtures.mk_part('P1', 500)
        sp1 = AstSlicePosition(250, PART_SLICE_POSTFIX_START, False)

        result = convert_slice_position_to_seq_position(part, sp1)

        self.assertEqual(result.x, 250)
        self.assertEqual(result.approximate, False)
        self.assertEqual(result.reference_sequence, reference)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEqual(ref, reference)
        self.assertEqual(abs_pos, 750)

    def test_convert_slice_position_to_seq_position_with_relative_to_end(self):

        reference, part = self.fixtures.mk_part('P2', 500)
        sp1 = AstSlicePosition(-125, PART_SLICE_POSTFIX_END, True)

        result = convert_slice_position_to_seq_position(part, sp1)

        self.assertEqual(result.x, -125)
        self.assertEqual(result.approximate, True)
        self.assertEqual(result.reference_sequence, reference)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEqual(ref, reference)
        self.assertEqual(abs_pos, 875)
