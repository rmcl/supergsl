from typing import List, Optional
from Bio.SeqRecord import SeqRecord
from supergsl.core.types import NucleotideSequence, PrimerPair


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

    def get_child_part_by_slice(self, identifier, start, end):
        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            start,
            end
        )
