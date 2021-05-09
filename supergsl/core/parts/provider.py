from typing import Optional
from supergsl.utils import import_class
from supergsl.core.exception import ConfigurationError
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition


class PartProvider(SuperGSLProvider):
    name : Optional[str] = None

    def __init__(self, name):
        self.name = name

    def get_provider_name(self):
        return self.name

    def list_parts(self):
        # The method is optional
        raise NotImplementedError(
            'List parts is not supported by "%s" part provider.' % self.name
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


class PartProviderPlugin(SuperGSLPlugin):
    name = 'part_provider'

    def register(self, symbol_table, compiler_settings):
        """Instantiate and register each part_providers defined in settings."""

        if 'part_providers' not in compiler_settings:
            raise ConfigurationError(
                'No part providers have been defined. Check your supergGSL settings.')

        import_symbol_table = symbol_table.nested_scope('imports')

        for provider_config in compiler_settings['part_providers']:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])
            provider_inst = provider_class(provider_config['name'], provider_config)

            provider_name = provider_inst.get_provider_name()
            if not provider_name:
                raise ConfigurationError('Provider "%s" does not specify a name.' % provider_class)

            import_symbol_table.insert(provider_name, provider_inst)
