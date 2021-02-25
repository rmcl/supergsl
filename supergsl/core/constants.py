"""A module of useful constants."""

THREE_PRIME = 'ThreePrime'
FIVE_PRIME = 'FivePrime'

PART_SLICE_POSTFIX_START = 'S'
PART_SLICE_POSTFIX_END = 'E'

UNAMBIGUOUS_DNA_SEQUENCE = 'DNA'
UNAMBIGUOUS_PROTEIN_SEQUENCE= 'PROTEIN'

# Import commonly used Sequence Onotology Terms
# Todo: Look into "tyto" (https://github.com/SynBioDex/tyto)
from sbol2.constants import (
    SO_GENE,
    SO_PROMOTER,
    SO_TERMINATOR,
    SO_CDS,
)

# SO_ORF seems not to be defined in pysbol3
SO_ORF = 'https://identifiers.org/SO:0000236'
SO_MRNA = 'https://identifiers.org/SO:0000234'
SO_HOMOLOGOUS_REGION = 'https://identifiers.org/SO:0000853'
SO_CDS_FRAGMENT = 'https://identifiers.org/SO:0001384'
