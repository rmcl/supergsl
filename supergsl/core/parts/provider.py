import re
from collections import namedtuple
from typing import Callable, Tuple
from supergsl.utils import import_class
from supergsl.core.exception import ConfigurationException
from supergsl.core.config import settings
from supergsl.core.exception import PartNotFoundException, ProviderNotFoundException
from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.parts import Part

class PartProvider(object):
    name = None

    def __init__(self, name):
        self.name = name

    def get_provider_name(self):
        return self.name

    def list_parts(self):
        # The method is optional
        raise NotImplementedError(
            'List parts is not supported by "%s" part provider.' % self.name
        )

    def resolve_import(self, identifier : str, alias : str) -> Tuple[re.Pattern, Callable[[str], Part]]:
        """Resolve the import of a part from this provider.

        Return a tuple with:
            * A regular expression to match symbols against
            * A callback method that given the actual identifier will return the `Part`.
        """
        pattern = re.compile(alias or identifier)
        return pattern, self.get_part

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part from the provider.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        raise NotImplementedError('Subclass to implement.')

    def get_child_part(self, identifier : str, sub_part_identifier : str) -> Part:
        raise NotImplementedError('Subclass to implement.')


class PartProviderPlugin(SuperGSLPlugin):
    name = 'part_provider'

    def register(self, symbol_table):
        self._initialize_part_providers(symbol_table)

    def _initialize_part_providers(self, symbol_table):
        self._providers = {}

        if 'part_providers' not in settings:
            raise ConfigurationException('No part providers have been defined. Check your supergGSL settings.')

        for provider_config in settings['part_providers']:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])
            provider_inst = provider_class(provider_config['name'], provider_config)

            provider_name = provider_inst.get_provider_name()
            if not provider_name:
                raise ConfigurationException('Provider "%s" does not specify a name.' % provider_class)

            symbol_table.register(provider_name, provider_inst)
