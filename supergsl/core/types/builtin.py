"""Define SuperGSL builtin types that can be used by plugins."""
from Bio.Seq import Seq
from typing import List, Type, Optional, Union
from collections import OrderedDict

from supergsl.core.sequence import SequenceEntry
from .base import SuperGSLType
from .position import Slice


class SuperGSLEnum(SuperGSLType):
    """Define a list of choices."""
    options : List[str] = []


class SliceInvertMixin:
    """A Mixin that allows a part to be sliced or inverted."""

    def slice(
        self,
        part_slice : Union[Slice, str],
        identifier : Optional[str] = None
    ) -> 'SuperGSLType':
        """Return a new part representing a sliced region."""
        raise Exception(f'{type(self)} does not support slicing.')

    def invert(self, identifier : Optional[str] = None):
        """Return a new Part with sequence on reverse strand."""
        raise Exception(f'{type(self)} does not support inverting.')


class NucleotideSequence(SuperGSLType, SliceInvertMixin):
    """A type representing a nucleotide sequence."""

    def __init__(self, sequence_entry : SequenceEntry):
        self.sequence_entry = sequence_entry

    @property
    def sequence(self) -> Seq:
        """Return the nucleotide sequence as a `Bio.Seq`."""
        return self.sequence_entry.sequence

    def __repr__(self):
        """Create a repr string representation of the Nucleotide Sequence."""
        return 'NucleotideSequence: {self.sequence}'


class AminoAcidSequence(SuperGSLType):
    """A type representing arbitrary an amino acid sequence."""

    def __init__(self, sequence_entry : SequenceEntry):
        self.sequence_entry = sequence_entry

    @property
    def sequence(self):
        """Return the amino acid sequence as a `Bio.Seq`."""
        return self.sequence_entry.sequence

    def __repr__(self):
        """Create a repr string representation of the Amino Acid Sequence."""
        return f'AminoAcidSequence: {self.sequence}'


class CodonTranslationTable(SuperGSLType):
    """Encode the mapping of DNA/RNA codons to protein.

    This is a thin wrapper around biopython's `Bio.Data.CodonTable`
    https://biopython.org/docs/1.75/api/Bio.Data.CodonTable.html
    """


class Collection(SuperGSLType):
    """Store a list of items."""

    def __init__(self, items : List[SuperGSLType]):
        self._items = items

    def __iter__(self):
        """Return an iterator over the items in the collection"""
        return iter(self._items)

    def __len__(self):
        """Return the number of items in the collection"""
        return len(self._items)

    def __repr__(self):
        """Create a repr string representation of the collection."""
        output = f'Collection (count: {len(self._items)})\n'
        for item_idx, item in enumerate(self._items):
            output += f'  {item_idx}. {item} \n'
        return output


class SliceAndInvertCollection(Collection):
    """A list of items that can be either inverted or sliced."""

    def __init__(
        self,
        item_collection : Collection,
        item_slice : Optional[Slice],
        inverted : bool
    ):
        super().__init__(item_collection._items)
        self.item_slice = item_slice
        self.inverted = inverted

    def __iter__(self):
        for item_symbol in self._items:
            symbol = item_symbol

            if self.item_slice:
                symbol = symbol.slice(self.item_slice)

            if self.inverted:
                symbol = symbol.invert()

            yield symbol
