import re
from re import Pattern, Match
from typing import Optional
from supergsl.utils import import_class
from supergsl.core.exception import ConfigurationError
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.parts import Part, SeqPosition

class PartProvider(object):
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

    def resolve_import(self, identifier : str, alias : str) -> Pattern:
        """Resolve the import of a part from this provider.

        Return a tuple with:
            * A regular expression to match symbols against
            * A callback method that given the actual identifier will return the `Part`.
        """
        return re.compile(alias or identifier)

    def get_symbol(self, identifier_match : Match):
        return self.get_part(identifier_match.string)

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part from the provider.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
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
        self.settings = compiler_settings
        self._initialize_part_providers(symbol_table)

    def _initialize_part_providers(self, symbol_table):
        self._providers = {}

        if 'part_providers' not in self.settings:
            raise ConfigurationException('No part providers have been defined. Check your supergGSL settings.')

        for provider_config in self.settings['part_providers']:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])
            provider_inst = provider_class(provider_config['name'], provider_config)

            provider_name = provider_inst.get_provider_name()
            if not provider_name:
                raise ConfigurationError('Provider "%s" does not specify a name.' % provider_class)

            symbol_table.register(provider_name, provider_inst)
