from typing import List, Optional
from Bio.SeqRecord import SeqRecord


class Part(object):
    """Represent a genomic piece of DNA."""

    def __init__(
        self,
        identifier,
        start_position,
        end_position,
        provider,
        parent_part=None,
        forward_primer=None,
        reverse_primer=None,
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

        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

        self.description = description
        self.alternative_names = alternative_names

        self._roles = []
        if roles:
            self._roles = roles


    @property
    def has_primers(self):
        return self.forward_primer and self.reverse_primer

    def set_primers(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

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

    def get_child_part_by_slice(self, identifier, start, end):
        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            start,
            end
        )
