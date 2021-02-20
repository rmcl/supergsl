from Bio.SeqRecord import SeqRecord
from supergsl.core.types import NucleotideSequence


class Part(NucleotideSequence):
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
        alternative_names=None
    ):
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

    @property
    def has_primers(self):
        return self.forward_primer and self.reverse_primer

    def set_primers(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer

    def get_sequence(self):
        ref, x = self.start.get_absolute_position_in_reference()
        ref_2, y = self.end.get_absolute_position_in_reference()

        if ref != ref_2:
            raise Exception("Reference sequences do not match.")

        description = self.description or ''
        return SeqRecord(
            ref[x:y],
            name=self.identifier,
            description=description)

    def get_child_part_by_slice(self, identifier, start, end):
        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            start,
            end
        )
