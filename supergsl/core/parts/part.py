from typing import List, Optional
from Bio.SeqRecord import SeqRecord
from supergsl.core.constants import (
    PART_SLICE_POSTFIX_START,
    PART_SLICE_POSTFIX_END
)
from supergsl.core.types import (
    SuperGSLType,
    NucleotideSequence,
    PrimerPair
)
from supergsl.core.ast import SlicePosition


class Part(NucleotideSequence):
    """Represent a genomic piece of DNA."""

    def __init__(
        self,
        identifier,
        start_position,
        end_position,
        provider,
        parent_part=None,
        extraction_primers=None,
        description=None,
        alternative_names=None,
        roles = None
    ):
        """Initialize a Part

        roles (optional): List of roles from standard ontology terms define the
            function of this part. Nice list of roles here:
                https://github.com/SynBioDex/pySBOL3/blob/master/sbol3/constants.py
        """
        self.identifier = identifier

        self.start = start_position
        self.end = end_position

        start_position.check_position_compatibility(end_position)

        self.provider = provider

        self.parent_part = parent_part

        self.extraction_primers = extraction_primers

        self.description = description
        self.alternative_names = alternative_names

        self._roles = []
        if roles:
            self._roles = roles

    @property
    def has_primers(self):
        """Returns true if this part has extraction primers defined."""
        return self.extraction_primers is not None

    def set_extraction_primers(self, extraction_primers : PrimerPair):
        """A PrimerPair which can be used to extract this Part from its template source."""
        self.extraction_primers = extraction_primers

    def get_extraction_primers(self) -> PrimerPair:
        """Return primers that will allow the extraction of this part from its template source."""
        return self.extraction_primers

    def get_sequence(self):
        ref, start_pos = self.start.get_absolute_position_in_reference()
        ref_2, end_pos = self.end.get_absolute_position_in_reference()

        if ref != ref_2:
            raise Exception("Reference sequences do not match.")

        description = self.description or ''
        return SeqRecord(
            ref[start_pos:end_pos],
            name=self.identifier,
            description=description)

    @property
    def roles(self):
        """Return a list of ontology terms for this part."""
        return self._roles

    def add_roles(self, roles : List[str]):
        """Add additional ontology terms to this part."""
        self._roles.extend(roles)
        self._roles = list(set(self._roles))

    def get_child_part_by_slice(self, identifier, start, end) -> 'Part':
        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            start,
            end
        )


    def convert_slice_position_to_seq_position(
        self,
        parent_part: 'Part',
        slice_position: SlicePosition
    ):
        """Convert SuperGSL `ast.SlicePosition` to `parts.SeqPosition` relative
        to the resolved parent part.

        ast.SlicePosition has the following properties:
            - index - The position on the part sequence.
            - postfix - "S" or "E" - Whether the index is relative to the start
                (S) or end (E) of the part.
            - approximate - Boolean flag determines if the index is approximate
                (True) or exact (False).
        """
        if slice_position.postfix == PART_SLICE_POSTFIX_START or slice_position.postfix is None:
            # the index is relative to the start of the part.
            return parent_part.start.get_relative_position(
                slice_position.index,
                slice_position.approximate)
        elif slice_position.postfix == PART_SLICE_POSTFIX_END:
            # the index is relative to the end of the part.
            return parent_part.end.get_relative_position(
                slice_position.index,
                slice_position.approximate)

        raise Exception('Unknown slice postfix: "%s"' % slice_position.postfix)


    def eval(self, ast_node: 'SymbolReference'):
        """Evaluate this part in the context of `ast.SymbolReference` node
        which is part of a SuperGSL Program"""

        if not ast_node.slice:
            # This node does not require slicing. no slice has been specified.
            return self

        if not ast_node.part:
            raise Exception('Node part not defined. Previous import pass failed.')

        parent_part = ast_node.part

        start = self.convert_slice_position_to_seq_position(parent_part, ast_node.slice.start)
        end = self.convert_slice_position_to_seq_position(parent_part, ast_node.slice.end)

        child_identifier = '%s[%s]' % (
            ast_node.identifier,
            ast_node.slice.get_slice_str()
        )
        ast_node.part = parent_part.get_child_part_by_slice(
            child_identifier, start, end)

        if ast_node.invert:
            raise NotImplementedError('Inverted parts not implemented yet!')

        return ast_node.part


class LazyLoadedPart(SuperGSLType):
    def eval(self, ast_node: 'SymbolReference') -> SuperGSLType:
        raise NotImplementedError('Subclass to implement.')
