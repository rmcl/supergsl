from Bio.Seq import Seq
from typing import List, Type
from collections import OrderedDict
from .base import SuperGSLType

class SuperGSLEnum(SuperGSLType):
    """Define a list of choices."""
    options : List[str] = []


class NucleotideSequence(SuperGSLType):
    """A type representing a nucleotide sequence."""

    def __init__(self, sequence : Seq):
        self._sequence = sequence

    def get_sequence(self) -> Seq:
        """Return the nucleotide sequence as a `Bio.Seq`."""
        return self._sequence

    def __repr__(self):
        return 'NucleotideSequence: %s' % self.get_sequence()


class AminoAcidSequence(SuperGSLType):
    """A type representing arbitrary an amino acid sequence."""

    def __init__(self, sequence : Seq):
        self._sequence = sequence

    def get_sequence(self):
        """Return the amino acid sequence as a `Bio.Seq`."""
        return self._sequence

    def __repr__(self):
        return 'AminoAcidSequence: %s' % self.get_sequence()

class CodonFrequencyTable(SuperGSLType):
    """
    https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/tree/master/codon_usage_data/tables
    """
    pass


class CodonTranslationTable(SuperGSLType):
    """Encode the mapping of DNA/RNA codons to protein.

    This is a thin wrapper around biopython's `Bio.Data.CodonTable`
    https://biopython.org/docs/1.75/api/Bio.Data.CodonTable.html
    """
    pass


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


class Collection(SuperGSLType):

    def __init__(self, items : List[SuperGSLType]):
        self._items = items

    def __iter__(self):
        return iter(self._items)

    def count(self):
        return len(self._items)

    def get_by_label(self, idx):
        pass

    def get_by_index(self, idx):
        pass

    def __repr__(self):
        output = 'Collection (count: %d)\n' % len(self._items)
        for item_idx, item in enumerate(self._items):
            output += '  %d. %s \n' % (item_idx, item)
        return output
