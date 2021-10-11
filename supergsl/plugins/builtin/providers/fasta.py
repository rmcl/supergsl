"""Implementation of a part provider backed by simple FASTA files."""
import csv
import gzip
from typing import Dict, List, Tuple, Any, TextIO
from mimetypes import guess_type
from Bio import SeqIO
from Bio.Seq import Seq

from supergsl.core.exception import PartNotFoundError
from supergsl.core.types.part import Part
from supergsl.core.types.slice import Slice
from supergsl.core.sequence import SequenceStore, SequenceEntry
from supergsl.core.parts import PartProvider, PartProviderConfig
from supergsl.plugins.pydna.primers import ExtractionPrimerBuilder

class FastaPartProvider(PartProvider):
    """Access parts provided by a simple FASTA file.

    Each entry in the FASTA file can be accessed as a part.

    Options:
    * fasta_file_path (required): The path to the FASTA file that will serve
        as part source
    * identifier_format (optional): A format string that will create identifiers
        from the record header of each fasta entry.
    """

    def __init__(self, name : str, config : PartProviderConfig):
        self._provider_name = name

        settings = config.provider_config
        self.fasta_file_path : str = settings['fasta_file_path']
        self.identifier_format : str = settings.get('identifier_format', '%s')

        self.sequence_store = config.sequence_store
        self._cached_parts : Dict[str, Part] = {}
        self._sequences_by_identifier : Dict[str, SequenceEntry] = {}
        self._loaded : bool = False

    def load(self) -> None:
        def _open(file_path : str) -> TextIO:
            encoding = guess_type(file_path)[1]
            if encoding == 'gzip':
                return gzip.open(file_path, mode='rt')
            else:
                return open(file_path, mode='rt')

        with _open(self.fasta_file_path) as fp:
            entries = SeqIO.parse(fp, 'fasta')

            self._sequences_by_identifier = {}
            for entry in entries:
                sequence_entry = self.sequence_store.add_from_reference(entry.seq)
                entry_name = self.identifier_format % entry.name
                self._sequences_by_identifier[entry_name] = sequence_entry

        self._loaded = True

    def list_parts(self):
        if not self._loaded:
            self.load()

        return list(self._sequences_by_identifier.keys())

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        try:
            return self._cached_parts[identifier]
        except KeyError:
            pass

        if not self._loaded:
            self.load()

        try:
            sequence_entry = self._sequences_by_identifier[identifier]
        except KeyError as error:
            raise PartNotFoundError('%s not found.' % identifier) from error

        part = Part(
            identifier,
            sequence_entry,
            provider=self)

        self._cached_parts[identifier] = part
        return part

    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        part_slice : Slice,
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""

        new_sequence_entry = self.sequence_store.slice(
            parent_part.sequence_entry,
            part_slice)

        return Part(
            identifier,
            sequence_entry=new_sequence_entry,
            provider=self,
            parent_part=parent_part
        )
