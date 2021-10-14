from typing import Dict, Any
from supergsl.core.constants import (
    FIVE_PRIME,
    THREE_PRIME
)
from .base import SuperGSLType

class AbsolutePosition:
    """Capture an absolute position in a sequence of fixed length."""

    def __init__(self, target_sequence_length, index, approximate):
        self.target_sequence_length = target_sequence_length
        self.index = index
        self.approximate = approximate

    @property
    def is_out_of_bounds(self) -> bool:
        """Return True if the index of this position exceeds the bounds of the sequence.

        For example, index=-15 or sequence_length + 25
        """
        if self.index < 0 or self.index >= self.target_sequence_length:
            return True
        return False

class Position:
    """Capture a position relative to a declared end of a Sequence."""

    def __init__(
        self,
        index: int,
        relative_to : str = FIVE_PRIME,
        approximate : bool = False
    ):
        self.index = index

        self.relative_to = relative_to

        # Position is relative to either the three prime or five prime side of the sequence.
        assert self.relative_to in [FIVE_PRIME, THREE_PRIME]

        self.approximate = approximate

    def compute_absolute_position(self, sequence_length : int) -> AbsolutePosition:
        """Given the length of the sequence compute the absolute position of this Position."""

        if self.relative_to == FIVE_PRIME:
            return AbsolutePosition(
                sequence_length,
                self.index,
                self.approximate)
        else:
            return AbsolutePosition(
                sequence_length,
                sequence_length - self.index - 1,
                self.approximate)


    def __repr__(self):
        return self.get_slice_pos_str()

    def serialize(self) -> Dict[str, Any]:
        """Serialize the `Position` object."""
        return {
            'index': self.index,
            'relative_to': self.relative_to,
            'approximate': self.approximate,
        }

    def get_slice_pos_str(self) -> str:
        """Return a string representation of the Position."""
        return '%s%d%s' % (
            '~' if self.approximate else '',
            self.index,
            '' if self.relative_to == FIVE_PRIME else 'E'
        )


class Slice(SuperGSLType):
    """Capture the start and end position of a part slice.

    If the end position is less than the start position then it is infered that
    the desired sequence lays on the reverse strand of a double-stranded molecule.
    """
    def __init__(self, start : Position, end : Position):
        self.start = start
        self.end = end


    @classmethod
    def from_five_prime_indexes(self, start_index, end_index):
        """Create a Slice from two sequence indexes both relative to the FIVE_PRIME side of the molecule."""
        start = Position(start_index)
        end = Position(end_index)

        return Slice(start, end)


    @classmethod
    def from_entire_sequence(self):
        """Create a Slice capturing the entire sequence."""
        start = Position(0, relative_to=FIVE_PRIME)
        end = Position(0, relative_to=THREE_PRIME)

        return Slice(start, end)


    @classmethod
    def from_str(cls, slice_string : str) -> 'Slice':
        """Create a `Slice` from a slice string.

        For example pGAL3[0:200S], the slice string would be "0:200S".
        """

        # Todo: Figure out how to resolve the circular dependency without putting
        # the import here.
        from supergsl.utils.slice import parse_slice_str
        return parse_slice_str(slice_string)

    def serialize(self) -> Dict[str, Any]:
        """Serialize the slice object."""
        return {
            'start': self.start.serialize(),
            'end': self.end.serialize()
        }

    def __repr__(self):
        return self.get_slice_str()

    def get_slice_str(self):
        """Return a string representation of the `Slice`."""
        return '%s:%s' % (
            self.start.get_slice_pos_str(),
            self.end.get_slice_pos_str()
        )
