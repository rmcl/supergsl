import typing
import re
from typing import Tuple, Callable, Optional, Mapping
from re import Pattern, Match
from supergsl.core.types import SuperGSLType
from supergsl.core.exception import PartSliceError
from supergsl.core.constants import (
    FIVE_PRIME,
    SO_GENE,
    SO_PROMOTER,
    SO_TERMINATOR,
    SO_HOMOLOGOUS_REGION,
    SO_CDS,
    SO_CDS_FRAGMENT,
    SO_MRNA
)

from supergsl.core.ast import SymbolReference
from supergsl.core.types.part import Part, LazyLoadedPart
from .provider import PartProvider


class UnknownPartPrefixError(Exception):
    pass


def get_promoter_len() -> int:
    """Get the configured length of a promoter region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM PROVIDER CONFIG
    """

    return 501

def get_terminator_len() -> int:
    return 501

def get_flank_len() -> int:
    """Get the configured length of a flanking region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM PROVIDER CONFIG
    """
    return 501

# Note: This is an unfortunate check here to make mypy feel good about mixins
# The parent should inherit from `supergsl.core.parts.PartProvider`, but
# that is rather difficult for mypy to infer when it runs.
if typing.TYPE_CHECKING:
    _Base = PartProvider
else:
    _Base = object


class PrefixedSliceLazyLoadedPart(LazyLoadedPart):
    def __init__(self, identifier, parent_part, prefix):
        self.identifier = identifier
        self.parent_part = parent_part
        self.part_prefix = prefix

    def eval(self) -> Part:
        return self.parent_part.get_prefixed_part(
            self.identifier,
            self.part_prefix)

    def __str__(self):
        'LazyPrefixPart: {} {}'.format(self.identifier, self.parent_part)

class PrefixedSlicePartProviderMixin(_Base):
    """A Part Mixin which enables support for fGSL prefixed-parts.

    From the GSL Paper, valid part prefixes are the following:

    g prefix gene locus gADH1 (equivalent to ORF prefix)
    p prefix promoter part pERG10
    t prefix terminator part tERG10
    u prefix upstream part uHO
    d prefix downstream part dHO
    o prefix open reading frame oERG10
    f prefix fusible ORF, no stop codon fERG10
    m prefix mRNA (ORF + terminator)

    Source:
        - https://github.com/Amyris/GslCore/blob/b738b3e107b91ed50a573b48d0dcf1be69c4ce6a/src/GslCore/CommonTypes.fs#L60
    """
    PART_TYPES = {
        'g': 'gene',
        'p': 'promoter',
        't': 'terminator',
        'u': 'upstream',
        'd': 'downstream',
        'o': 'orf',
        'f': 'fusible_orf',
        'm': 'mRNA'
    }

    def resolve_import(
        self,
        identifier : str,
        alias : Optional[str]
    ) -> Mapping[str, SuperGSLType]:
        """Resolve the import of a part from this provider.

        Override the default functionality and instead register a part for each
        of the part prefixes available to the user.
        """

        new_symbols = {}
        part_identifier = alias or identifier
        for part_prefix in self.PART_TYPES.keys():
            part_name = '{}{}'.format(part_prefix, part_identifier)
            lazy_part = PrefixedSliceLazyLoadedPart(
                identifier,
                self,
                part_prefix)

            new_symbols[part_name] = lazy_part

        return new_symbols


    def get_prefixed_part(self, identifier : str, prefix : str) -> Part:
        parent_part = self.get_part(identifier)
        if prefix == '':
            return parent_part

        try:
            part_type = self.PART_TYPES[prefix]
        except KeyError:
            raise UnknownPartPrefixError('Invalid part prefix "{}" for "{}".'.format(
                prefix,
                identifier,

            ))

        start_pos, end_pos = self.build_part_type_slice_pos(parent_part, part_type)
        roles = self.get_roles_by_part_type(part_type)
        part = self.get_child_part_by_slice(
            parent_part=parent_part,
            identifier='{}{}'.format(prefix, identifier),
            start=start_pos,
            end=end_pos)
        part.add_roles(roles)
        return part


    def get_roles_by_part_type(self, part_type):
        """Return a list of roles based on part type."""
        part_type_role_map = {
            'gene': [SO_GENE],
            'promoter': [SO_PROMOTER],
            'terminator': [SO_TERMINATOR],
            'upstream': [SO_HOMOLOGOUS_REGION],
            'downstream': [SO_HOMOLOGOUS_REGION],
            'orf': [SO_CDS],
            'fusible_orf': [SO_CDS_FRAGMENT],
            'mRNA': [SO_MRNA],
        }

        return part_type_role_map[part_type]

    def build_part_type_slice_pos(self, parent_part, part_slice_type):
        """Build the slice of a part based on the requested part type.

        parts often have a part type specified by the prefix, for example
        p for promoter.

        Refer to translateGenePrefix in GslCore for reference logic
        https://github.com/Amyris/GslCore/blob/d2c613907d33b110a2f53021146342234e0d8f3b/src/GslCore/DnaCreation.fs#L53

        """
        if part_slice_type == 'promoter':
            new_start = parent_part.start.get_relative_position(
                x=-1 * get_promoter_len(),
                approximate=True)

            return new_start, parent_part.start

        elif part_slice_type == 'orf':
            return parent_part.start, parent_part.end

        elif part_slice_type == 'gene':
            return parent_part.start, parent_part.end

        elif part_slice_type == 'upstream':
            new_start = parent_part.start.get_relative_position(
                x=-1 * get_flank_len(),
                approximate=True)

            return new_start, parent_part.start

        elif part_slice_type == 'downstream':
            new_end = parent_part.end.get_relative_position(
                x=get_flank_len(),
                approximate=True)

            return parent_part.end, new_end

        elif part_slice_type == 'terminator':
            new_end = parent_part.end.get_relative_position(
                x=get_terminator_len(),
                approximate=True)

            return parent_part.end, new_end

        raise PartSliceError('"%s" prefix is not implemented yet.' % part_slice_type)
