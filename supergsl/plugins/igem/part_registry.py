"""Retrieve BioBrick parts."""
from typing import cast
from Bio.Seq import Seq
from Bio import SeqIO
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.role import SO_HOMOLOGOUS_REGION
from supergsl.core.parts.provider import ConstantPartProvider

from supergsl.plugins.builtin.providers.synbiohub import SynBioHubPartProvider
from .types import BioBrickPart
from .constants import (
    BIOBRICK_PREFIX_SEQUENCE,
    BIOBRICK_SUFFIX_SEQUENCE
)

class PartRegistry(SynBioHubPartProvider):
    """Retrieve BioBrick parts from the iGEM Part Registry.

    Subclass the SynBioHubPartProvider and create BioBrick parts.
    """

    def get_part(self, identifier : str) -> BioBrickPart:
        """Retrieve a BioBrick part by synbiohub identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `BioBrickPart`
        """

        part_details = self.get_part_details(identifier)
        part = BioBrickPart.from_payload_sequence(
            sequence_store=self.sequence_store,
            payload_sequence=part_details['sequence'],
            identifier=identifier,
            provider=self,
            description=part_details['description'],
            roles=part_details['roles'])

        return part


class BioBrickPartProvider(ConstantPartProvider):
    """A Part provider for returning constant regions from the BioBrick standard.

    More about the BioBrick standard can be found here:
        http://parts.igem.org/Help:Standards/Assembly/RFC10

    To configure this part provider add the following to `supergsl-config.json`:
    ```
    {
        "name": "biobrick",
        "provider_class": "supergsl.plugins.igem.BioBrickPartProvider"
    }

    References
    * http://parts.igem.org/Help:Standards/Assembly/RFC10
    """

    PART_DETAILS = {
        'prefix': (
            'Biobrick RFC[10] prefix',
            BIOBRICK_PREFIX_SEQUENCE,
            [SO_HOMOLOGOUS_REGION]
        ),
        'suffix': (
            'Biobrick RFC[10] suffix',
            BIOBRICK_SUFFIX_SEQUENCE,
            [SO_HOMOLOGOUS_REGION]
        )
    }
