from Bio.Seq import Seq
from supergsl.core.sequence import SequenceEntry

from .base import SuperGSLType
from .builtin import NucleotideSequence


class Primer(NucleotideSequence):
    """Represent a short nucleotide sequence that provides a starting point for replication"""
    pass


class PairedPrimer(Primer):
    """A primer with a `body` and `tail` region useful for annealing two parts."""
    def __init__(self, body : Seq, tail : Seq):
        self._body = body
        self._tail = tail

    @property
    def sequence(self) -> Seq:
        return Seq(str(self._tail) + str(self._body))

    @property
    def body(self) -> Seq:
        """Return a Sequence representing "body" region of the primer."""
        return self._body

    @property
    def tail(self) -> Seq:
        """Return a Sequence representing "tail" region of the primer."""
        return self._tail


class PrimerPair(SuperGSLType):
    """A pair of primers used in PCR reactions to amplify a region of DNA."""

    @classmethod
    def from_sequences(cls, forward_primer_seq, reverse_primer_seq):
        """Construct a PrimerPair from two sequences."""
        return PrimerPair(Primer(forward_primer_seq), Primer(reverse_primer_seq))

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
