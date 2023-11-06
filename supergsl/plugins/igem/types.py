from Bio.Seq import Seq
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.position import Position

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

    @classmethod
    def from_payload_sequence(
        cls,
        sequence_store,
        payload_sequence,
        **kwargs
    ) -> 'BioBrickPart':

        check_is_valid_biobrick(payload_sequence)
        sequence = Seq(''.join([
            str(BIOBRICK_PREFIX_SEQUENCE).lower(),
            str(payload_sequence).lower(),
            str(BIOBRICK_SUFFIX_SEQUENCE).lower()
        ]))

        sequence_entry = sequence_store.add_from_reference(sequence)

        return BioBrickPart(
            sequence_entry=sequence_entry,
            **kwargs
        )
