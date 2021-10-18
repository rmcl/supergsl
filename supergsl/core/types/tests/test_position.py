from unittest import TestCase
from supergsl.core.constants import FIVE_PRIME, THREE_PRIME
from supergsl.core.types.position import (
    AbsoluteSlice,
    AbsolutePosition,
    Position,
    Slice
)


class AbsolutePositionTestCase(TestCase):
    """Test the AbsolutePosition class."""

    def test_derive_from_relative_five_prime_position(self):
        """Compute a new absolute position using a relative position defined relative to FIVE PRIME."""
        start = AbsolutePosition(1000, 100, False)
        rel_pos = Position(50, FIVE_PRIME, False)

        new_abs_pos = start.derive_from_relative_position(rel_pos)
        self.assertEqual(new_abs_pos.target_sequence_length, 1000)
        self.assertEqual(new_abs_pos.index, 150)

    def test_derive_from_relative_three_prime_position(self):
        """Compute a new absolute position using a relative position defined relative to THREE PRIME."""
        start = AbsolutePosition(1000, 100, False)
        rel_pos = Position(50, THREE_PRIME, False)

        new_abs_pos = start.derive_from_relative_position(rel_pos)
        self.assertEqual(new_abs_pos.target_sequence_length, 1000)
        self.assertEqual(new_abs_pos.index, 150)


class AbsoluteSliceTestCase(TestCase):
    """Test the AbsoluteSlice class."""

    def test_derive_from_relative_slice_five_prime(self):
        """Get the correct absolute position when adjusted for a relative position."""
        abs_slice = AbsoluteSlice(
            AbsolutePosition(1000, 100, False),
            AbsolutePosition(1000, 400, False))

        new_abs_slice = abs_slice.derive_from_relative_slice(
            Slice.from_five_prime_indexes(50,100)
        )

        self.assertEqual(new_abs_slice.start.index, 150)
        self.assertEqual(new_abs_slice.end.index, 200)

    def test_derive_from_relative_slice_reverse_strand(self):
        abs_slice = AbsoluteSlice(
            AbsolutePosition(1000, 100, False),
            AbsolutePosition(1000, 400, False))

        new_abs_slice = abs_slice.derive_from_relative_slice(
            Slice.from_five_prime_indexes(100,50)
        )

        self.assertEqual(new_abs_slice.start.index, 200)
        self.assertEqual(new_abs_slice.end.index, 150)

    def test_derive_from_relative_slice_three_prime(self):
        abs_slice = AbsoluteSlice(
            AbsolutePosition(1000, 100, False),
            AbsolutePosition(1000, 400, False))

        new_abs_slice = abs_slice.derive_from_relative_slice(Slice(
            Position(50, FIVE_PRIME),
            Position(-50, THREE_PRIME))
        )

        self.assertEqual(new_abs_slice.start.index, 150)
        self.assertEqual(new_abs_slice.end.index, 350)
