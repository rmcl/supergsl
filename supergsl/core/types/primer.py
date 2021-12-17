from Bio.Seq import Seq
from supergsl.core.sequence import SequenceEntry, Role
from supergsl.utils.sequence import filter_links_by_roles
from .base import SuperGSLType
from .builtin import NucleotideSequence


paired_primer_body_role = Role(
    'supergsl/builtin/primer/body',
    'Primer Body',
    (
        'The region of a paired primer which has homology to the particular '
        'sequence or part to be amplified.'
    )
)
paired_primer_tail_role = Role(
    'supergsl/builtin/primer/tail',
    'Primer Tail',
    'The region of a paired primer which hangs off of a part.'
)

class Primer(NucleotideSequence):
    """Represent a short nucleotide sequence that provides a starting point for replication"""
    pass


class PairedPrimer(Primer):
    """A primer with a `body` and `tail` region useful for annealing two parts."""
    def __init__(self, primer_sequence_entry : SequenceEntry):
        self.primer_sequence_entry = primer_sequence_entry

    @property
    def sequence(self) -> Seq:
        return self.primer_sequence_entry.sequence

    @property
    def body(self) -> Seq:
        """Return a Sequence representing "body" region of the primer."""
        links = filter_links_by_roles(self.primer_sequence_entry, [paired_primer_body_role])
        if len(links) == 0:
            raise Exception('Paired Primer does not have a body sub-sequence.')
        if len(links) > 1:
            raise Exception('Paired Primer has more than one body sub-sequence.')

        return links[0].sequence

    @property
    def tail(self) -> Seq:
        """Return a Sequence representing "tail" region of the primer."""

        links = filter_links_by_roles(self.primer_sequence_entry, [paired_primer_tail_role])
        if len(links) == 0:
            raise Exception('Paired Primer does not have a tail sub-sequence.')
        if len(links) > 1:
            raise Exception('Paired Primer has more than one tail sub-sequence.')

        return links[0].sequence


class PrimerPair(SuperGSLType):
    """A pair of primers used in PCR reactions to amplify a region of DNA."""

    @classmethod
    def from_sequence_entries(
        cls,
        forward_primer_sequence : SequenceEntry,
        reverse_primer_sequence : SequenceEntry
    ):
        """Construct a PrimerPair from two sequences."""
        return PrimerPair(Primer(forward_primer_sequence), Primer(reverse_primer_sequence))

    def __init__(self, forward_primer, reverse_primer):
        self._forward_primer : Primer = forward_primer
        self._reverse_primer : Primer = reverse_primer

    @property
    def forward(self) -> Primer:
        """Return the forward primer of the pair"""
        return self._forward_primer

    @property
    def reverse(self) -> Primer:
        """Return the reverse primer of the pair"""
        return self._reverse_primer
