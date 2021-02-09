from re import Pattern
from typing import Optional, Tuple, List
from collections import OrderedDict
from supergsl.core.exception import SymbolNotFoundException
from supergsl.core.exception import ProviderNotFoundException

class SymbolTable(object):

    def __init__(self):
        self._path_to_plugins = {}
        self._symbols_providers: List[Tuple[Pattern, object]] = []
        self._symbols = {}

    def register(self, plugin_import_path: str, plugin_class: object) -> None:
        """Register a plugin provider to be available for import a given import path."""
        self._path_to_plugins[plugin_import_path] = plugin_class

    def get_plugin_provider(self, plugin_import_path):
        provider = self._path_to_plugins.get(plugin_import_path, None)
        if not provider:
            raise ProviderNotFoundException('Provider for %s not found.' % plugin_import_path)

        return provider

    def resolve_symbol(
        self,
        plugin_import_path: str,
        import_name: str,
        alias: Optional[str] = None
    ):
        """Resolve a symbol from an import statement of a superGSL program."""

        print('Attempting to resolve symbol: %s, %s, %s' % (plugin_import_path, import_name, alias))

        active_alias = alias or import_name
        try:
            self.get_symbol(active_alias)
        except SymbolNotFoundException:
            pass
        else:
            raise Exception('Alias "%s" imported twice.' % active_alias)

        plugin_provider = self.get_plugin_provider(plugin_import_path)
        alias_pattern = plugin_provider.resolve_import(import_name, active_alias)

        self._symbols_providers.append((alias_pattern, plugin_provider))
        return self.get_symbol(active_alias)

    def get_symbol(self, identifier):
        """Retrieve an imported symbol."""

        try:
            return self._symbols[identifier]
        except KeyError:
            pass

        for symbol_pattern, symbol_handler in self._symbols_providers:
            match = symbol_pattern.fullmatch(identifier)
            if match:
                self._symbols[identifier] = symbol_handler.get_symbol(match)
                return self._symbols[identifier]

        raise SymbolNotFoundException('"%s" is not defined.' % identifier)
