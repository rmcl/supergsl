"""Retrieve Parts from a MixedLibrary."""
from typing import cast, Dict, List, Optional, Union
from Bio.Seq import Seq
from Bio import SeqIO
from supergsl.core.types import SuperGSLType
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.builtin import Collection
from supergsl.core.types.protein import Protein
from supergsl.core.provider import ProviderConfig
from supergsl.core.parts.provider import PartProvider

from mixed_parts import open_library
from mixed_parts.exceptions import RoleDoesNotExist

class MixedPartLibraryProvider(PartProvider):
    """Retrieve and Save Parts to the Mixed Part Library File Format.

    """

    def __init__(self, name : str, config : ProviderConfig):
        self._provider_name = name
        self._cached_parts: Dict[str, Part] = {}
        self.sequence_store = config.sequence_store
        self._config = config

        settings = config.settings

        self._library = open_library(settings['library_path'])

    @property
    def library(self):
        return self._library

    def _create_part_from_details(self, part_details : dict) -> Union[Part, Protein]:
        """Create a SuperGSL part from details provided."""
        sequence_entry = self.sequence_store.add_from_reference(
            part_details['sequence'])

        if part_details['part_type'] == 'protein':
            return Protein(
                identifier=part_details['identifier'],
                sequence_entry=sequence_entry,
                provider=self,
                description=part_details['description'],
                roles=part_details['roles'])
        elif part_details['part_type'] == 'dna':
            return Part(
                identifier=part_details['identifier'],
                sequence_entry=sequence_entry,
                provider=self,
                description=part_details['description'],
                roles=part_details['roles'])

        raise ValueError('Unknown part type "%s"' % part_details['part_type'])



    def create_collection(self, name : str, description : str) -> Collection:
        pass

    def add_to_collection(self, collection : Collection, part : Part):
        pass

    def get_collection(self, name : str) -> Collection:
        collection_detail = self._library.collections.get(name)

        collection_items = []
        for item in collection_detail['items']:
            item_part = self.get(item['part_id'])
            collection_items.append(item_part)

        return Collection(collection_items)

    def list_collections(self) -> List[str]:
        return [
            collection_detail['name']
            for collection_detail in self._library.collections.list()
        ]


    def get(self, identifier : str) -> SuperGSLType:
        """Retrieve a Part from mixed part library by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Union[Part]`
        """

        part_details = self._library.parts.get(identifier)
        return self._create_part_from_details(part_details)


    def save(self, part : Part):
        """Save a part to the mixed part library"""

        for role in part.roles:
            try:
                self._library.parts.get_role(role.uri)
            except RoleDoesNotExist:
                self._library.parts.create_role(
                    uri=role.uri,
                    name=role.name,
                    description=role.description
                )

        self._library.parts.create_part(
            identifier=part.identifier,
            part_type='dna',
            sequence=str(part.sequence),
            description=part.description,
            roles=[
                role.uri
                for role in part.roles
            ]
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
        library_part_details = self._library.parts.list(
            part_type='dna',
            include_sequence=True)
        parts = []

        for part_details in library_part_details:
            # Todo: Filter
            part = self._create_part_from_details(part_details)
            parts.append(part)

        return parts
