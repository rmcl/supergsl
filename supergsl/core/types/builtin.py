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


class CodonTranslationTable(SuperGSLType):
    """Encode the mapping of DNA/RNA codons to protein.

    This is a thin wrapper around biopython's `Bio.Data.CodonTable`
    https://biopython.org/docs/1.75/api/Bio.Data.CodonTable.html
    """
    pass


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
