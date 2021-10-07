"""Implement access to SynBioHub powered genetic part repos."""
from typing import Dict
from Bio.Seq import Seq
from sbol2 import Document, PartShop
from supergsl.core.exception import ConfigurationError
from supergsl.core.constants import THREE_PRIME
from supergsl.utils.cache import FileCache

from supergsl.core.sequence import SequenceStore
from supergsl.core.parts import PartProvider, PartProviderConfig
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition


class SynBioHubPartProvider(PartProvider):
    """A Part provider for accessing SynBioHub powered genetic part repos.

    From the pysbol2 docs here here are two examples:
        https://synbiohub.org/public/igem
        https://synbiohub.org

    To configure this part provider add the following to `supergsl-config.json`:
    ```
    {
        "name": "igem",
        "provider_class": "supergsl.plugins.builtin.providers.SynBioHubPartProvider",
        "repository_url": "https://synbiohub.org/public/igem",
        "repository_username": null,
        "repository_password": null
    }
    ```
    `repository_username` and `repository_password` are optional.

    References
    * https://pysbol.readthedocs.io/en/latest/repositories.html
    * https://synbiohub.org/public/igem
    """

    def __init__(self, name : str, config : PartProviderConfig):
        self._provider_name = name
        self._cached_parts: Dict[str, Part] = {}

        settings = config.provider_config
        self.repository_url = settings.get('repository_url', None)
        if not self.repository_url:
            ConfigurationError('"%s" requires that repository_url be set.')

        self.enable_part_cache = settings.get('enable_part_cache', True)
        self._cache = FileCache(name, enable = self.enable_part_cache)

        self.repository_username = settings.get('repository_username', None)
        self.repository_password = settings.get('repository_password', None)

    def get_part_details(self, identifier : str) -> dict:
        """Retrieve Part details from the remote repository."""

        try:
            return self._cache.get(identifier)
        except KeyError:
            pass

        part_doc = Document()
        part_shop = PartShop(self.repository_url)
        if self.repository_username and self.repository_password:
            part_shop.login(self.repository_username, self.repository_password)

        part_shop.pull(identifier, part_doc)

        component_definition = part_doc.componentDefinitions[identifier]

        details = {
            'roles': component_definition.roles,
            'description': component_definition.description,
            'sequence': Seq(component_definition.compile())
        }

        self._cache.store(identifier, details)
        return details

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

        part_details = self.get_part_details(identifier)
        reference_sequence = part_details['sequence']

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
            provider=self,
            description=part_details['description'],
            roles=part_details['roles'])

        self._cached_parts[identifier] = part
        return part
