import unittest
import mock
from Bio.Seq import Seq
from supergsl.core.ast import (
    Slice as AstSlice,
    SlicePosition as AstSlicePosition,
    Part as AstPart
)
from supergsl.core.parts import Part, SeqPosition

from supergsl.core.test.fixtures import SuperGSLCoreFixtures
from supergsl.core.parts.slice import ResolvePartSlicePass
from supergsl.core.constants import (
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)


class ResolvePartSlicePassTestCase(unittest.TestCase):
    """Test case for `ResolvePartSlicePass`."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_convert_slice_position_to_seq_position_with_relative_to_start(self):
        #index: int, postfix : str, approximate : bool

        reference, part = self.fixtures.mk_part('P1', 500)
        sp1 = AstSlicePosition(250, PART_SLICE_POSTFIX_START, False)

        rps = ResolvePartSlicePass(None)
        result = rps.convert_slice_position_to_seq_position(part, sp1)

        self.assertEquals(result.x, 250)
        self.assertEquals(result.approximate, False)
        self.assertEquals(result.parent, part.start)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEquals(ref, reference)
        self.assertEquals(abs_pos, 750)

    def test_convert_slice_position_to_seq_position_with_relative_to_end(self):
        #index: int, postfix : str, approximate : bool

        reference, part = self.fixtures.mk_part('P2', 500)
        sp1 = AstSlicePosition(-125, PART_SLICE_POSTFIX_END, True)

        rps = ResolvePartSlicePass(None)
        result = rps.convert_slice_position_to_seq_position(part, sp1)

        self.assertEquals(result.x, -125)
        self.assertEquals(result.approximate, True)
        self.assertEquals(result.parent, part.end)

        ref, abs_pos = result.get_absolute_position_in_reference()
        self.assertEquals(ref, reference)
        self.assertEquals(abs_pos, 875)

