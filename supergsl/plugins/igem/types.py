from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition

from .constants import (
    BIOBRICK_PREFIX_SEQUENCE,
    BIOBRICK_SUFFIX_SEQUENCE
)
from .utils import check_is_valid_biobrick


class BioBrickPart(Part):
    """A Part consistent with the BioBrick Standard.

    From the docs, "BioBrick RFC[10] it must not contain the following
    restriction sites, as these are unique to the prefix and suffix"
    Taken from http://parts.igem.org/Help:Standards/Assembly/RFC10.
    """
    def __init__(self, payload_sequence, **kwargs):
        check_is_valid_biobrick(payload_sequence)
        self._payload_sequence = payload_sequence

        sequence = ''.join([
            str(BIOBRICK_PREFIX_SEQUENCE).lower(),
            str(payload_sequence).lower(),
            str(BIOBRICK_SUFFIX_SEQUENCE).lower()
        ])

        kwargs['start_position'] = SeqPosition.from_reference(
            x=0,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=sequence
        )
        kwargs['end_position'] = kwargs['start_position'].get_relative_position(
            x=len(sequence))

        super().__init__(**kwargs)
