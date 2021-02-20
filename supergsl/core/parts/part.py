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
        alternative_names=None
    ):
        self.identifier = identifier

        self.start = start_position
        self.end = end_position

        start_position.check_position_compatibility(end_position)

        self.provider = provider

        self.parent_part = parent_part

        self.extraction_primers = extraction_primers

        self.description = description
        self.alternative_names = alternative_names

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

    def get_child_part_by_slice(self, identifier, start, end):
        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            start,
            end
        )
