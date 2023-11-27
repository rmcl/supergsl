"""Retrieve BioBrick parts."""
from typing import cast, Dict, List, Optional
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

    def _create_part_from_details(self, part_details : dict) -> Part:
        """Create a SuperGSL part from details provided."""
        sequence_entry = self.sequence_store.add_from_reference(
            part_details['sequence'])

        part = Part(
            identifier=part_details['identifier'],
            sequence_entry=sequence_entry,
            provider=self,
            description=part_details['description'],
            roles=part_details['roles'])

        return part

    def get_part(self, identifier : str) -> Part:
        """Retrieve a Part from mixed part library by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        part_details = self._library.get(identifier)
        return self._create_part_from_details(part_details)



    def save_part(self, part : Part):
        """Save a part to the mixed part library"""
        self._library.create_part(
            identifier=part.identifier,
            part_type='dna',
            sequence=str(part.sequence),
            description=part.description,
            roles=part.roles
        )

    def search(
        self,
        query : Optional[str] = None,
        roles : Optional[List[str]] = None
    ) -> List[Part]:
        """Search for parts in the part library.

        Arguments:
            query: match identifier or description
            roles: list[str] match parts with a given role
        """
        library_part_details = self._library.list(
            part_type='dna',
            include_sequence=True)
        parts = []

        for part_details in library_part_details:
            # Todo: Filter
            part = self._create_part_from_details(part_details)
            parts.append(part)

        return parts
