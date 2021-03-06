"""Define SuperGSL builtin types that can be used by plugins."""
from Bio.Seq import Seq
from typing import List


class SuperGSLType(object):
    """Base class defining types available in SuperGSL."""
    pass


class NucleotideSequence(SuperGSLType):
    """A type representing a nucleotide sequence."""

    def get_sequence(self) -> Seq:
        """Return the nucleotide sequence as a `Bio.Seq`."""
        raise NotImplementedError('Subclass to implement.')


class AminoAcidSequence(SuperGSLType):
    """A type representing arbitrary an amino acid sequence."""

    def get_sequence(self) -> Seq:
        """Return the amino acid sequence as a `Bio.Seq`."""
        raise NotImplementedError('Subclass to implement.')


class Primer(NucleotideSequence):
    """Represent a short nucleotide sequence that provides a starting point for replication"""
    def __init__(self, primer_seq : Seq):
        self._sequence = primer_seq

    def get_sequence(self) -> Seq:
        """Return the amino acid sequence as a `Bio.Seq`."""
        return self._sequence


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
