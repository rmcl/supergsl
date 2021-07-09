import csv
import gzip
from typing import Dict, List, Tuple, Any, TextIO
from mimetypes import guess_type
from Bio import SeqIO
from Bio.Seq import Seq

from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.position import SeqPosition
from supergsl.core.exception import PartNotFoundError
from supergsl.core.types.part import Part
from supergsl.core.parts import PartProvider
from supergsl.core.parts.prefix_part import PrefixedSlicePartProviderMixin
from supergsl.plugins.pydna.primers import ExtractionPrimerBuilder

class FastaPartProvider(PartProvider):
    """Access parts provided by a simple FASTA file."""

    def __init__(self, name : str, settings : dict):
        self._provider_name = name
        self.fasta_file_path : str = settings['fasta_file_path']
        self.identifier_format : str = settings.get('identifier_format', '%s')
        self._cached_parts : Dict[str, Part] = {}
        self._sequences_by_entry : Dict[str, Seq] = {}
        self._loaded = False

    def load(self) -> None:

        def _open(file_path : str) -> TextIO:
            encoding = guess_type(file_path)[1]
            if encoding == 'gzip':
                return gzip.open(file_path, mode='rt')
            else:
                return open(file_path, mode='rt')

        with _open(self.fasta_file_path) as fp:
            entries = SeqIO.parse(fp, 'fasta')
            print(entries)

            self._sequences_by_entry = {}
            for entry in entries:
                entry_name = self.identifier_format % entry.name
                self._sequences_by_entry[entry_name] = entry

        self._loaded = True

    def list_parts(self):
        if not hasattr(self, '_sequences_by_entry'):
            self.load()

        return list(self._sequences_by_entry.keys())

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
            reference_sequence = self._sequences_by_entry[identifier].seq
        except KeyError as error:
            raise PartNotFoundError('%s not found.' % identifier) from error

        start = SeqPosition.from_reference(
            x=0,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=reference_sequence
        )

        end = start.get_relative_position(
            x=len(reference_sequence))

        part = Part(
            identifier,
            start,
            end,
            provider=self)

        self._cached_parts[identifier] = part
        return part

    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        start : SeqPosition,
        end : SeqPosition
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""

        return Part(
            identifier,
            start,
            end,
            provider=self,
            parent_part=parent_part
        )
