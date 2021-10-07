from typing import Dict, Any
from supergsl.core.constants import (
    FIVE_PRIME,
    THREE_PRIME,
    STRAND_WATSON,
    STRAND_CRICK
)
from .base import SuperGSLType


class Position:
    """Capture a position.

    This differs from `supergsl.core.types.position.SeqPosition` in that it does
    not store the associated reference sequence."""

    def __init__(
        self,
        index: int,
        relative_to : str = FIVE_PRIME,
        approximate : bool = False,
        strand : str = STRAND_WATSON
    ):
        self.index = index

        self.relative_to = relative_to

        # Position is relative to either the three prime or five prime side of the sequence.
        assert self.relative_to in [FIVE_PRIME, THREE_PRIME]

        self.approximate = approximate
        self.strand = strand

        # Watson strand or crick
        # TODO: DEFINE THESE AS CONSTANTS in constants.py
        assert self.strand in [STRAND_WATSON, STRAND_CRICK]

    def __repr__(self):
        return self.get_slice_pos_str()

    def serialize(self) -> Dict[str, Any]:
        """Serialize the `Position` object."""
        return {
            'index': self.index,
            'relative_to': self.relative_to,
            'approximate': self.approximate,
            'strand': self.strand
        }

    def get_slice_pos_str(self) -> str:
        """Return a string representation of the Position."""
        return '%s%d%s' % (
            '~' if self.approximate else '',
            self.index,
            '' if self.relative_to == FIVE_PRIME else 'E'
        )


class Slice(SuperGSLType):
    """Capture the start and end position of a part slice."""
    def __init__(self, start : Position, end : Position):
        self.start = start
        self.end = end

        assert self.start.strand == self.end.strand, 'Strands of start and end positions must match'

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
