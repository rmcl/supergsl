from supergsl.core.ast import SlicePosition
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition
from supergsl.core.constants import (
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)


def convert_slice_position_to_seq_position(
    parent_part : Part,
    slice_position : SlicePosition
) -> SeqPosition:
    """Convert SuperGSL `ast.SlicePosition` to `types.position.SeqPosition` relative
    to the parent part.

    ast.SlicePosition has the following properties:
        - index - The position on the part sequence.
        - postfix - "S" or "E" - Whether the index is relative to the start
            (S) or end (E) of the part.
        - approximate - Boolean flag determines if the index is approximate
            (True) or exact (False).
    """
    if slice_position.postfix == PART_SLICE_POSTFIX_START or slice_position.postfix is None:
        # the index is relative to the start of the part.
        return parent_part.start.get_relative_position(
            slice_position.index,
            slice_position.approximate)

    elif slice_position.postfix == PART_SLICE_POSTFIX_END:
        # the index is relative to the end of the part.
        return parent_part.end.get_relative_position(
            slice_position.index,
            slice_position.approximate)

    raise Exception('Unknown slice postfix: "%s"' % slice_position.postfix)
