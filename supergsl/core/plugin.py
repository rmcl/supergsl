"""Support for SuperGSL's plugin infrastructure."""

import inspect
import importlib
from typing import Dict
from supergsl.core.config import settings
from supergsl.core.exception import ConfigurationException
from supergsl.core.symbol_table import SymbolTable

class SuperGSLPlugin(object):
    """Base class for defining a SuperGSL Plugin."""

    def register(self, symbol_table):
        """Register Functions, enums, etc that the plugin provides.

        Example: symbol_table.register(import_path, mod_class)
        """
        pass


class PluginProvider(object):
    """Resolve and register SuperGSL plugins."""

    def __init__(self, symbol_table: SymbolTable):
        self._plugins: Dict[str, SuperGSLPlugin] = {}
        self._symbol_table: SymbolTable = symbol_table

        if 'plugins' not in settings:
            raise ConfigurationException(
                'No plugins have been defined. Check your supergGSL settings.')

        for plugin_path in settings['plugins']:
            print('Resolving plugin "%s"' % plugin_path)

            self.resolve_plugins_from_config(plugin_path)


    def resolve_plugins_from_config(self, module_path: str) -> None:
        """Attempt to resolve and register a plugin at a specific path."""
        module = importlib.import_module(module_path)
        module_classes = inspect.getmembers(module, inspect.isclass)

        for _, plugin_class in module_classes:
            if issubclass(plugin_class, SuperGSLPlugin):
                print('Registering plugin...', plugin_class)

                plugin_inst = plugin_class()
                plugin_inst.register(self._symbol_table)
