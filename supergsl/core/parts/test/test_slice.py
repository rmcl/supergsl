import unittest
import mock
from Bio.Seq import Seq
from supergsl.core.ast import (
    Slice as AstSlice,
    SlicePosition as AstSlicePosition,
    SymbolReference as AstSymbolReference
)
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition

from supergsl.core.test.fixtures import SuperGSLCoreFixtures
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
        self.assertEqual(result.parent, part.start)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEqual(ref, reference)
        self.assertEqual(abs_pos, 750)

    def test_convert_slice_position_to_seq_position_with_relative_to_end(self):

        reference, part = self.fixtures.mk_part('P2', 500)
        sp1 = AstSlicePosition(-125, PART_SLICE_POSTFIX_END, True)

        result = convert_slice_position_to_seq_position(part, sp1)

        self.assertEqual(result.x, -125)
        self.assertEqual(result.approximate, True)
        self.assertEqual(result.parent, part.end)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEqual(ref, reference)
        self.assertEqual(abs_pos, 875)

    """
    def test_visit_part_node(self):
        _, pre_slice_part = self.fixtures.mk_part('P3', 500)
        pre_slice_part.provider.get_child_part_by_slice.return_value = 'SLICED-PART'

        start = AstSlicePosition(-125, PART_SLICE_POSTFIX_END, True)
        end = AstSlicePosition(0, PART_SLICE_POSTFIX_END, True)
        slice_node = AstSlice(start, end)

        part_node = AstSymbolReference('TEST', slice_node, False)
        part_node.part = pre_slice_part

        rps = ResolvePartSlicePass(None)
        new_part_node = rps.visit_part_node(part_node)

        pre_slice_part.provider.get_child_part_by_slice.assert_called_once()
        create_child_args = pre_slice_part.provider.get_child_part_by_slice.call_args.args

        self.assertEquals(create_child_args[0], pre_slice_part)
        self.assertEquals(create_child_args[1], 'TEST[~-125E:~0E]')

        part_start = create_child_args[2]
        self.assertEquals(part_start.x, -125)
        self.assertEquals(part_start.approximate, True)

        part_end = create_child_args[3]
        self.assertEquals(part_end.x, 0)
        self.assertEquals(part_end.approximate, True)

        self.assertEquals(new_part_node.part, 'SLICED-PART')
        self.assertEquals(part_node, new_part_node)
        self.assertEquals(new_part_node.parent_parts, [pre_slice_part])
    """
