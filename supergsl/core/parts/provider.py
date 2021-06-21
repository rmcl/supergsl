from typing import Optional, Dict
from Bio.Seq import Seq

from supergsl.utils import import_class

from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition

from supergsl.core.exception import ConfigurationError, PartNotFoundError
from supergsl.core.constants import THREE_PRIME


class PartProvider(SuperGSLProvider):

    def __init__(self, name, compiler_settings):
        self._provider_name = name
        self._compiler_settings = compiler_settings

    @property
    def provider_name(self):
        """Return the name of this part provider."""
        return self._provider_name

    def list_parts(self):
        """Return all parts available through this provider."""
        raise NotImplementedError(
            'List parts is not supported by "%s" part provider.' % self.provider_name
        )

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Resolve a part from the provider and register it in the symbol table."""
        part_identifier = alias or identifier
        symbol_table.insert(part_identifier, self.get_part(identifier))

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part from the provider.

        Arguments:
            identifier  A identifier to select a part from this provider
        """
        raise NotImplementedError('Subclass to implement.')

    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        start : SeqPosition,
        end : SeqPosition
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""
        raise NotImplementedError('Subclass to implement')


class ConstantPartProvider(PartProvider):
    """Base class for providing constant parts.

    To utilize this class, subclass and define a PART_DETAILS dictionary.
    ```
    {
        <part-identifier>: (
            <description>
            <sequence>
            [<list of roles>]
        )
    }
    ```
    """

    PART_DETAILS = {
        'part-name': ('part description', 'ATGC', ['LIST OF ROLES']),
    }

    def __init__(self, name : str, compiler_settings : dict):
        self._provider_name = name
        self._settings = compiler_settings
        self._cached_parts: Dict[str, Part] = {}


    def get_part_details(self, part_identifier):
        """Return constant details about a part."""
        return self.PART_DETAILS[part_identifier]

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

        try:
            description, reference_sequence, roles = \
                self.get_part_details(identifier)

        except KeyError as error:
            raise PartNotFoundError(
                'Part "%s" not found.' % identifier) from error

        if not isinstance(reference_sequence, Seq):
            reference_sequence = Seq(reference_sequence)

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
            description=description,
            roles=roles)

        self._cached_parts[identifier] = part
        return part


class PartProviderPlugin(SuperGSLPlugin):
    """Plugin stub to register part providers."""
    name = 'part_provider'

    def register(self, compiler_settings):
        """Instantiate and register each part_providers defined in settings."""

        if 'part_providers' not in compiler_settings:
            raise ConfigurationError(
                'No part providers have been defined. Check your supergGSL settings.')

        for provider_config in compiler_settings['part_providers']:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])
            provider_inst = provider_class(provider_config['name'], provider_config)

            self.register_provider(provider_config['name'], provider_inst)
