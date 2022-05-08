import typing
import re
from typing import Tuple, Callable, Optional, Mapping, Dict
from re import Pattern, Match
from supergsl.core.types import SuperGSLType
from supergsl.core.types.position import Slice, Position
from supergsl.core.exception import PartSliceError
from supergsl.core.constants import (
    FIVE_PRIME,
    THREE_PRIME
)

from supergsl.core.types.role import (
    GENE,
    PROMOTER,
    TERMINATOR,
    HOMOLOGOUS_REGION,
    CDS,
    CDS_FRAGMENT,
    MRNA
)

from supergsl.core.types.part import Part
from .provider import PartProvider


class UnknownPartPrefixError(Exception):
    pass


def get_promoter_len() -> int:
    """Get the configured length of a promoter region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM PROVIDER CONFIG
    """

    return 501

def get_terminator_len() -> int:
    """Get the configured length of a terminator region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM PROVIDER CONFIG
    """
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

        # Todo: these are implemented in fGSL, but more research is needed to
        # understand what their implementation looks like.
        #'f': 'fusible_orf',
        #'m': 'mRNA'
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

        new_symbols : Dict[str, SuperGSLType] = {}
        part_identifier = alias or identifier

        # Add the parent part
        parent_part = self.get_part(identifier)
        new_symbols[part_identifier] = parent_part

        # Add all the prefix parts.
        for part_prefix in self.PART_TYPES.keys():
            prefixed_part_name = f'{part_prefix}{part_identifier}'
            new_symbols[prefixed_part_name] = self.get_prefixed_part(
                identifier,
                part_prefix)

        return new_symbols


    def get_prefixed_part(self, identifier : str, prefix : str) -> Part:
        parent_part = self.get_part(identifier)
        if prefix == '':
            return parent_part

        try:
            part_type = self.PART_TYPES[prefix]
        except KeyError as error:
            raise UnknownPartPrefixError(
                f'Invalid part prefix "{prefix}" for "{identifier}"') from error

        new_slice = self.build_part_type_slice_pos(parent_part, part_type)

        roles = self.get_roles_by_part_type(part_type)
        part = self.get_child_part_by_slice(
            parent_part=parent_part,
            identifier=f'{prefix}{identifier}',
            part_slice=new_slice)
        part.add_roles(roles)
        return part


    def get_roles_by_part_type(self, part_type):
        """Return a list of roles based on part type."""
        part_type_role_map = {
            'gene': [GENE],
            'promoter': [PROMOTER],
            'terminator': [TERMINATOR],
            'upstream': [HOMOLOGOUS_REGION],
            'downstream': [HOMOLOGOUS_REGION],
            'orf': [CDS],
            'fusible_orf': [CDS_FRAGMENT],
            'mRNA': [MRNA],
        }

        return part_type_role_map[part_type]

    def build_part_type_slice_pos(self, parent_part : Part, part_slice_type : str) -> Slice:
        """Build the slice of a part based on the requested part type.

        parts often have a part type specified by the prefix, for example
        p for promoter.

        Refer to translateGenePrefix in GslCore for reference logic
        https://github.com/Amyris/GslCore/blob/d2c613907d33b110a2f53021146342234e0d8f3b/src/GslCore/DnaCreation.fs#L53

        """
        if part_slice_type == 'promoter':
            return Slice(
                Position(-1 * get_promoter_len(), relative_to=FIVE_PRIME, approximate=True),
                Position(0, relative_to=FIVE_PRIME))

        if part_slice_type in ['orf', 'gene']:
            return Slice.from_entire_sequence()

        if part_slice_type == 'upstream':
            return Slice(
                Position(-1 * get_flank_len(), relative_to=FIVE_PRIME, approximate=True),
                Position(0, relative_to=FIVE_PRIME))

        if part_slice_type == 'downstream':
            return Slice(
                Position(0, relative_to=THREE_PRIME),
                Position(1 * get_flank_len(), relative_to=THREE_PRIME, approximate=True))

        if part_slice_type == 'terminator':
            return Slice(
                Position(0, relative_to=THREE_PRIME),
                Position(get_terminator_len(), relative_to=THREE_PRIME, approximate=True))

        raise PartSliceError(f'"{part_slice_type}" prefix is not implemented yet.')
