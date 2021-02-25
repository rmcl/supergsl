"""Define a hierarchy of types that can be used by SuperGSL and its plugins."""
from .builtin import (
    SuperGSLType,
    SuperGSLEnum,

    NucleotideSequence,
    AminoAcidSequence,

    CodonFrequencyTable,
    CodonTranslationTable
)

TYPES = [
    SuperGSLType,
    SuperGSLEnum,
    NucleotideSequence,
    AminoAcidSequence,
    CodonFrequencyTable,
    CodonTranslationTable
]

def resolve_type(type_name):
    """Look up a SuperGSL type by name."""
    for type in TYPES:
        if type.__class__ == type_name:
            return type

    raise Exception('Unknown Type "%s"' % type_name)
