from typing import List, Optional, Dict, Union
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from supergsl.core.exception import PartError
from supergsl.core.sequence import SequenceEntry
from supergsl.core.types import SuperGSLType
from supergsl.core.types.builtin import NucleotideSequence
from supergsl.core.types.primer import PrimerPair
from supergsl.core.types.position import Slice


class Part(NucleotideSequence):
    """Represent a genomic piece of DNA."""

    def __init__(
        self,
        identifier : str,
        sequence_entry : SequenceEntry,
        provider,
        parent_part : Optional['Part'] = None,
        extraction_primers : Optional[PrimerPair] = None,
        description : Optional[str] = None,
        alternative_names : Optional[List[str]] = None,
        roles : Optional[List[str]]= None
    ):
        """Initialize a Part

        roles (optional): List of roles from standard ontology terms define the
            function of this part. Nice list of roles here:
                https://github.com/SynBioDex/pySBOL3/blob/master/sbol3/constants.py
        """
        self.identifier = identifier

        self.sequence_entry = sequence_entry

        self.provider = provider

        self.parent_part = parent_part

        self.extraction_primers = extraction_primers

        self.description = description
        self.alternative_names = alternative_names

        self._roles = []
        if roles:
            self._roles = roles

    @property
    def has_primers(self):
        """Returns true if this part has extraction primers defined."""
        return self.extraction_primers is not None

    def set_extraction_primers(self, extraction_primers : PrimerPair):
        """A PrimerPair which can be used to extract this Part from its template source."""
        self.extraction_primers = extraction_primers

    def get_extraction_primers(self) -> PrimerPair:
        """Return primers that will allow the extraction of this part from its template source."""
        if not self.extraction_primers:
            raise PartError('Extraction Primers for part have not been defined.')

        return self.extraction_primers

    @property
    def sequence(self) -> Seq:
        return self.sequence_entry.sequence

    @property
    def sequence_record(self) -> SeqRecord:
        """Return a `Bio.SeqRecord` entry for this part."""
        description = self.description or ''
        return SeqRecord(
            self.sequence,
            name=self.identifier,
            description=description)

    @property
    def roles(self):
        """Return a list of ontology terms for this part."""
        return self._roles

    def add_roles(self, roles : List[str]):
        """Add additional ontology terms to this part."""
        self._roles.extend(roles)
        self._roles = list(set(self._roles))

    def slice(self, part_slice : Union[Slice, str], identifier : Optional[str] = None) -> 'Part':
        """Return a new part representing a sliced region.

        part_slice is either a Slice object or a str.
        """
        if isinstance(part_slice, str):
            part_slice = Slice.from_str(part_slice)

        if not identifier:
            identifier = '%s[%s]' % (
                self.identifier,
                part_slice.get_slice_str()
            )

        return self.provider.get_child_part_by_slice(
            self,
            identifier,
            part_slice
        )

    def eval(self):
        """Evaluate this part."""
        return self

    def __repr__(self):
        return self.identifier

    def serialize(self) -> Dict:
        return {
            'identifier': self.identifier,
            'description': self.description,
            'sequence': str(self.sequence),
        }

    def print(self) -> str:
        """Display details about the SuperGSL object."""
        return '%s: %s' % (self.identifier, self.sequence)


class LazyLoadedPart(SuperGSLType):
    def eval(self) -> SuperGSLType:
        raise NotImplementedError('Subclass to implement.')
