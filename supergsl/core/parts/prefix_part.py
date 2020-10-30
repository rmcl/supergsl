from supergsl.core.constants import FIVE_PRIME
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
    return 500

def get_flank_len() -> int:
    """Get the configured length of a flanking region.

    ### TODO: REFACTOR THIS METHOD TO COME FROM PROVIDER CONFIG
    """
    return 500

class PrefixedPartSliceMixin(object):
    """A Part Mixin which enables support for fGSL prefixed-parts."""

    def get_child_part(self, part_prefix, alias=None):
        try:
            part_type = self.get_part_type(part_prefix)
        except UnknownPartPrefixError:
            return super().get_child_part(part_prefix)

        start_pos, end_pos = self.build_part_type_slice_pos(part_type)
        print('PREFIX', start_pos, end_pos)
        child_part = self.get_child_part_by_slice(
            parent_part=self,
            identifier=(alias or part_prefix),
            start=start_pos,
            end=end_pos)

        return child_part


    def get_part_type(self, prefix) -> str:
        """Validate the part prefix.

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

        try:
            return PART_TYPES[prefix]
        except KeyError:
            raise UnknownPartPrefixError('Invalid part prefix "%s" in "%s".' % (
                prefix,
                self.identifier
            ))

    def build_part_type_slice_pos(self, part_type):
        """Build the slice of a part based on the requested part type.

        parts often have a part type specified by the prefix, for example
        p for promoter.

        Refer to translateGenePrefix in GslCore for reference logic
        https://github.com/Amyris/GslCore/blob/d2c613907d33b110a2f53021146342234e0d8f3b/src/GslCore/DnaCreation.fs#L53

        """
        if part_type == 'promoter':
            new_start = self.start.get_relative_position(
                x=-1 * get_promoter_len(),
                approximate=True)

            return new_start, self.start

        elif part_type == 'gene':
            return self.start, self.end

        elif part_type == 'upstream':
            new_start = self.start.get_relative_position(
                x=-1 * get_flank_len(),
                approximate=True)

            return new_start, self.start

        elif part_type == 'downstream':
            new_end = self.end.get_relative_position(
                x=get_flank_len(),
                approximate=True)

            return self.end, new_end

        elif part_type == 'terminator':
            new_end = self.end.get_relative_position(
                x=get_terminator_len(),
                approximate=True)

            return self.end, new_end

        raise Exception('"%s" not implemented yet.' % part_type)

class PrefixedPart(PrefixedPartSliceMixin, Part):
    pass
