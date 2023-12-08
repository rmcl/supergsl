from typing import List, Optional
from supergsl.core.sequence import SequenceEntry
from supergsl.core.types.builtin import AminoAcidSequence
from supergsl.utils.sequence import build_truncated_sequence

class Protein(AminoAcidSequence):
    """A type representing arbitrary an amino acid sequence."""

    def __init__(
        self,
        identifier : str,
        sequence_entry : SequenceEntry,
        alternative_names : Optional[List[str]] = None,
        description : Optional[str] = None,
        roles : Optional[List[str]] = None
    ):
        """Initialize a Protein."""
        self.identifier = identifier
        self.sequence_entry = sequence_entry
        self.alternative_names = alternative_names
        self.description = description
        self.roles = roles

    def __repr__(self):
        """Create a repr string representation of the Amino Acid Sequence."""
        return f'Protein {self.identifier}: {build_truncated_sequence(self.sequence)}'
