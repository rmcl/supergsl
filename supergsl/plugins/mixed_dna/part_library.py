"""Retrieve BioBrick parts."""
from typing import cast, Dict
from Bio.Seq import Seq
from Bio import SeqIO
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.provider import ProviderConfig
from supergsl.core.parts.provider import PartProvider

from mixed_parts import open_library


class MixedPartLibraryProvider(PartProvider):
    """Retrieve and Save Parts to the Mixed Part Library File Format.

    """

    def __init__(self, name : str, config : ProviderConfig):
        self._provider_name = name
        self._cached_parts: Dict[str, Part] = {}
        self.sequence_store = config.sequence_store

        settings = config.settings

        self._library = open_library(settings['library_path'])

    def get_part(self, identifier : str) -> Part:
        """Retrieve a Part from mixed part library by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        part_details = self._library(identifier)
        sequence_entry = self.sequence_store.add_from_reference(
            part_details['sequence'])

        part = Part(
            identifier=identifier,
            sequence_entry=sequence_entry,
            provider=self,
            description=part_details['description'],
            roles=part_details['roles'])

        return part

    def save_part(self, part : Part):
        """Save a part to the mixed part library"""
        self._library.create_part(
            identifier=part.identifier,
            part_type='dna',
            sequence=str(part.sequence),
            description=part.description,
            roles=part.roles
        )
