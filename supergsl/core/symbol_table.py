from re import Pattern
from typing import Optional
from collections import OrderedDict
from supergsl.core.exception import SymbolNotFoundException

class SymbolTable(object):

    def __init__(self):
        self._path_to_plugins = {}
        self._symbols : Tuple[Pattern, ] = []

    def register(self, plugin_import_path : str, plugin_class) -> None:
        """Register a plugin provider to be available for import a given import path."""
        self._path_to_plugins[plugin_import_path] = plugin_class

    def get_plugin_provider(self, plugin_import_path):
        print(self._path_to_plugins)
        provider = self._path_to_plugins.get(plugin_import_path, None)
        if not provider:
            raise Exception('No plugin PROVIDER!', plugin_import_path)

        return provider

    def resolve_symbol(
        self,
        plugin_import_path : str,
        import_name : str,
        alias : Optional[str] = None
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

        self._symbols.append((alias_pattern, plugin_provider))
        return self.get_symbol(active_alias)

    def get_symbol(self, identifier):
        """Retrieve an imported symbol."""

        for symbol_pattern, symbol_handler in self._symbols:
            match = symbol_pattern.fullmatch(identifier)
            if match:
                return symbol_handler.get_symbol(match)

        raise SymbolNotFoundException('"%s" is not defined.' % identifier)
