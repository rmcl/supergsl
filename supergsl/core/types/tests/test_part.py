from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.constants import (
    SO_PROMOTER,
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)
from supergsl.core.ast import SlicePosition as AstSlicePosition
from supergsl.core.types.part import convert_slice_position_to_seq_position



class PartTestCase(TestCase):
    """Test case for SeqPosition."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_part_roles_are_recapitulated(self):
        _, p1 = self.fixtures.mk_part('SWEET', 300, roles=[
            SO_PROMOTER
        ])

        self.assertEqual(p1.roles, [SO_PROMOTER])

    def test_part_print(self):
        """Test the part print statement outputs in a nice format."""
        _, part = self.fixtures.mk_part(identifier='BOOM', part_seq_len=25)

        self.assertEqual(
            part.print(),
            "BOOM: %s" % part.sequence)

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
