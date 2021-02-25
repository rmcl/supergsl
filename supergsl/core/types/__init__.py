from .builtin import (
    SuperGSLType,

    NucleotideSequence,
    AminoAcidSequence,

    Primer,
    PrimerPair,
)

TYPES = [
    SuperGSLType,
    NucleotideSequence,
    AminoAcidSequence,

    Primer,
    PrimerPair,
]

def resolve_type(type_name):
    for type in TYPES:
        if type.__class__ == type_name:
            return type

    raise Exception('Unknown Type "%s"' % type_name)
