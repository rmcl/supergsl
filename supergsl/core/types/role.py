from typing import NamedTuple


class Role(NamedTuple):
    """Represent a sequence role."""
    uri : str
    name : str
    description : str

    def __eq__(self, other):
        """Roles with the same uri should be considered identical."""
        return self.uri == other.uri


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

GENE = Role(SO_GENE, 'Gene', '')
PROMOTER = Role(SO_PROMOTER, 'Promoter', '')
TERMINATOR = Role(SO_TERMINATOR, 'Terminator', '')

ORF = Role(SO_ORF, 'Open Reading Frame', '')
MRNA = Role(SO_MRNA, 'Messenger RNA', '')
HOMOLOGOUS_REGION = Role(SO_HOMOLOGOUS_REGION, 'Homologous Region', '')
CDS = Role(SO_CDS, 'Coding DNA Sequence', '')
CDS_FRAGMENT = Role(SO_CDS_FRAGMENT, 'Coding DNA Sequence Fragment', '')

role_to_biopython_types = {
    GENE: 'gene',
    PROMOTER: 'promoter',
    TERMINATOR: 'terminator',
    ORF: 'orf',
    CDS_FRAGMENT: 'cds'
}

def convert_role_to_biopython_type(role : Role) -> str:
    """Convert a given `Role` to the corresponding biopython type."""
    try:
        return role_to_biopython_types[role]
    except TypeError as error:
        raise Exception(
            'Unknown mapping from {role.uri} to biopython type.') from error
