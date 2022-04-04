from Bio.Seq import Seq
from supergsl.core.sequence import SequenceStore, SequenceEntry, Role, SequenceAnnotation
from supergsl.utils.sequence import filter_annotations_by_roles
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

    @classmethod
    def from_sequences(self, sequence_store : SequenceStore, body : Seq, tail : Seq):
        """Create a PairedPrimer from body and tail sequences."""
        primer_seq_entry = sequence_store.add_from_reference(
            body + tail,
            annotations=[
                SequenceAnnotation.from_five_prime_indexes(
                    start=0,
                    end=len(body),
                    roles=[paired_primer_body_role],
                    payload={}),
                SequenceAnnotation.from_five_prime_indexes(
                    start=len(body),
                    end=len(body) + len(tail),
                    roles=[paired_primer_tail_role],
                    payload={})
            ])

        return PairedPrimer(primer_seq_entry)

    @property
    def sequence(self) -> Seq:
        return self.primer_sequence_entry.sequence

    @property
    def body(self) -> Seq:
        """Return a Sequence representing "body" region of the primer."""

        # TODO: SHOULD WE FIGURE OUT HOW TO CACHE SLICING THE SAME PART MULTIPLE TIMES!?!?!?

        annotations = filter_annotations_by_roles(
            self.primer_sequence_entry,
            [paired_primer_body_role])
        if len(annotations) == 0:
            raise Exception('Paired Primer does not have a body sub-sequence.')
        if len(annotations) > 1:
            raise Exception('Paired Primer has more than one body sub-sequence.')

        sequence_store = self.primer_sequence_entry.sequence_store
        body_entry = sequence_store.slice_from_annotation(
            self.primer_sequence_entry,
            annotations[0])

        return body_entry.sequence

    @property
    def tail(self) -> Seq:
        """Return a Sequence representing "tail" region of the primer."""

        annotations = filter_annotations_by_roles(
            self.primer_sequence_entry,
            [paired_primer_tail_role])
        if len(annotations) == 0:
            raise Exception('Paired Primer does not have a tail sub-sequence.')
        if len(annotations) > 1:
            raise Exception('Paired Primer has more than one tail sub-sequence.')

        sequence_store = self.primer_sequence_entry.sequence_store
        body_entry = sequence_store.slice_from_annotation(
            self.primer_sequence_entry,
            annotations[0])

        return body_entry.sequence


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
