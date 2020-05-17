from supergsl.utils import import_class
from supergsl.core.exception import ConfigurationException
from supergsl.core.config import settings
from supergsl.core.backend import BreadthFirstNodeFilteredPass


class PartSymbolTable(object):

    def __init__(self):
        self._parts = {}

        self._initialize_providers()

    def get_part(self, part_alias):
        try:
            return self._parts[part_alias]
        except KeyError:
            raise Exception('Part "%s" has not been defined.' % part_alias)

    def resolve_part(self, provider_name, part_name, alias=None):
        print('Resolving Part: %s, %s, %s' % (provider_name, part_name, alias))
        if not alias:
            alias = part_name

        if alias in self._parts:
            raise Exception('Part "%s" is already defined.' % alias)

        provider = self.resolve_provider(provider_name)
        part = provider.get_part(part_name)

        self._parts[alias] = part

    def _initialize_providers(self):
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

            self._providers[provider_name] = provider_inst

    def resolve_provider(self, provider_name):
        try:
            return self._providers[provider_name]
        except KeyError:
            raise Exception('Unknown part provider "%s".' % provider_name)


class Part(object):
    def __init__(self, name, sequence):
        self.name = None
        self.sequence = None


class PartProvider(object):
    name = None

    def __init__(self, name):
        self.name = name

    def get_provider_name(self):
        return self.name

    def get_part(self, identifier):
        """Retrieve a part from the provider.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """
        raise NotImplemented('Subclass to implement.')


class ResolvePartPass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'ProgramImport': self.visit_import_node,
            'Part': self.visit_part_node,
        }

    def before_pass(self, ast):
        self.part_symbol_table = PartSymbolTable()
        ast.symbol_registry.register('parts', self.part_symbol_table)
        return ast

    def visit_import_node(self, node):
        for program_import in node.imports:
            self.part_symbol_table.resolve_part(
                node.module,
                program_import.identifier
            )

    def visit_part_node(self, node):
        part = self.part_symbol_table.get_part(node.part_name)
        node.part = part
