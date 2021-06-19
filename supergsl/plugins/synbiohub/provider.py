from typing import Dict
from Bio.Seq import Seq
from sbol2 import Document, PartShop
from supergsl.core.exception import ConfigurationError
from supergsl.core.constants import THREE_PRIME

from supergsl.core.parts import PartProvider
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
        "provider_class": "supergsl.plugins.synbiohub.SynBioHubPartProvider",
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

    def __init__(self, name : str, settings : dict):
        self.name = name
        self._cached_parts: Dict[str, Part] = {}
        self.repository_url = settings.get('repository_url', None)
        if not self.repository_url:
            ConfigurationError('"%s" requires that repository_url be set.')

        self.repository_username = settings.get('repository_username', None)
        self.repository_password = settings.get('repository_password', None)

    def retrieve_part_details(self, identifier : str) -> dict:
        """Retrieve Part details from the remote repository."""
        part_doc = Document()
        part_shop = PartShop(self.repository_url)
        if self.repository_username and self.repository_password:
            part_shop.login(self.repository_username, self.repository_password)

        part_shop.pull(identifier, part_doc)

        component_definition = part_doc.componentDefinitions[identifier]

        return {
            'roles': component_definition.roles,
            'description': component_definition.description,
            'sequence': Seq(component_definition.compile())
        }

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

        part_details = self.retrieve_part_details(identifier)
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
