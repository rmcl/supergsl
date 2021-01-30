import typing
from typing import Tuple, Callable
import re
from re import Pattern, Match
from supergsl.core.constants import FIVE_PRIME
from supergsl.core.exception import PartSliceError
from .part import Part
from .position import SeqPosition

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


class PrefixedSlicePartProviderMixin(object):
    """A Part Mixin which enables support for fGSL prefixed-parts."""

    # Note: This is an unfortunate check here to make mypy feel good about using
    # self.identifier which is typically set in the parent class which this is
    # mixed into.
    if typing.TYPE_CHECKING:
        identifier = None

    def resolve_import(self, identifier : str, alias : str) -> Tuple[Pattern, Callable[[str], Part]]:
        """Resolve the import of a part from this provider."""

        pattern = self.get_prefix_pattern(alias or identifier)

        def get_part_handler(identifier : str, pattern_match : Match):
            matched_identifier = pattern_match.group('identifier')
            matched_prefix = pattern_match.group('prefix')

            parent_part = self.get_part(matched_identifier)
            if matched_prefix == '':
                return parent_part
            else:
                part_type = self.get_part_type(matched_prefix)
                start_pos, end_pos = self.build_part_type_slice_pos(parent_part, part_type)
                return self.get_child_part_by_slice(
                    parent_part=parent_part,
                    identifier=identifier,
                    start=start_pos,
                    end=end_pos)

            print(identifier, matched_identifier, matched_prefix)

            return self.get_child_part(
                matched_prefix,
                alias=identifier)

        return pattern, get_part_handler

    def get_prefix_pattern(self, identifier):
        allowed_prefixes = ''.join(self.PART_TYPES.keys())
        return re.compile('(?P<prefix>[%s]?)(?P<identifier>%s)' % (
            allowed_prefixes,
            identifier
        ))

    """
    https://github.com/Amyris/GslCore/blob/b738b3e107b91ed50a573b48d0dcf1be69c4ce6a/src/GslCore/CommonTypes.fs#L60

    From the GSL Paper, valid part prefixes are the following:
    g prefix gene locus gADH1 (equivalent to ORF prefix)
    p prefix promoter part pERG10
    t prefix terminator part tERG10
    u prefix upstream part uHO
    d prefix downstream part dHO
    o prefix open reading frame oERG10
    f prefix fusible ORF, no stop codon fERG10
    m prefix mRNA (ORF + terminator)
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

    def get_part_type(self, prefix) -> str:
        """Validate the part prefix."""
        try:
            return self.PART_TYPES[prefix]
        except KeyError:
            raise UnknownPartPrefixError('Invalid part prefix "%s" in "%s".' % (
                prefix,
                self.identifier
            ))

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

        raise PartSliceError('"%s" prefix is not implemented yet.' % part_type)
