"""Support for SuperGSL's plugin infrastructure."""

import inspect
import importlib
from typing import Dict, Set, Type, Mapping, cast
from supergsl.core.exception import (
    ConfigurationError,
    NotFoundError,
    SymbolNotFoundError
)
from supergsl.core.types import SuperGSLType
from supergsl.core.provider import SuperGSLProvider, ProviderGroup
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.sequence import SequenceStore


class SuperGSLPlugin(object):
    """Base class for defining a SuperGSL Plugin.

    An example of writing a function for SuperGSL:

    class ExamplePlugin(SuperGSLPlugin):

        def register(self, compiler_settings):
            self.register_function(
                'import_path',
                'functionname',
                SuperGSLFunctionDeclaration(FunctionClass, settings))
    """

    def __init__(self, symbol_table : SymbolTable, compiler_settings : dict):
        self.symbol_table = symbol_table
        self.functions : Dict[str, SuperGSLFunctionDeclaration] = {}
        self.register(compiler_settings)

    def resolve_import(
        self,
        identifier : str,
        alias : str
    ) -> Mapping[str, SuperGSLType]:
        """Import a identifier and register it in the symbol table."""
        if identifier not in self.functions:
            raise NotFoundError('%s not found in module.' % identifier)

        symbol_alias = alias or identifier
        return {
            symbol_alias: self.functions[identifier]
        }

    def get_or_create_provider_group_for_module_path(self, module_path):
        import_symbol_table = self.symbol_table.enter_nested_scope('imports')
        try:
            provider_group = import_symbol_table.lookup(module_path)
        except SymbolNotFoundError:
            provider_group = ProviderGroup()
            import_symbol_table.insert(module_path, provider_group)

        return provider_group

    def register_provider(
        self,
        import_path : str,
        provider_inst : SuperGSLProvider
    ):
        """Register a provider allowing the user to import things from a SuperGSLProgram."""
        provider_group = self.get_or_create_provider_group_for_module_path(import_path)
        provider_group.add_provider(provider_inst)

    def register_function(
        self,
        import_path : str,
        function_name : str,
        function_declaration : SuperGSLFunctionDeclaration
    ):
        """Register a function making it available for import in SuperGSL."""

        # Todo: Improve the setup of function declarations so this is better.
        function_declaration.set_sequence_store(
            cast(SequenceStore, self.symbol_table.lookup('sequences')))

        self.functions[function_name] = function_declaration

        provider_group = self.get_or_create_provider_group_for_module_path(import_path)
        if self not in provider_group:
            provider_group.add_provider(self)


    def register(self, compiler_settings : dict):
        """Register Functions, enums, etc that the plugin provides.

        Example: symbol_table.register(import_path, mod_class)
        """


class PluginProvider(object):
    """Resolve and register SuperGSL plugins."""

    def __init__(self, symbol_table: SymbolTable, compiler_settings : dict):
        # We add SuperGSLPlugin base type so that we don't attempt to register it.
        self._resolved_plugins : Set[Type] = set([SuperGSLPlugin])
        self._symbol_table: SymbolTable = symbol_table
        self._compiler_settings = compiler_settings

        if 'plugins' not in compiler_settings:
            raise ConfigurationError(
                'No plugins have been defined. Check your supergGSL settings.')

        for plugin_path in compiler_settings['plugins']:
            print('Resolving plugin path "%s"' % plugin_path)

            self.resolve_plugins_from_config(plugin_path)

    def get_adhoc_plugin(self):
        return SuperGSLPlugin(self._symbol_table, self._compiler_settings)


    def resolve_plugins_from_config(self, module_path: str) -> None:
        """Attempt to resolve and register a plugin at a specific path."""

        try:
            module = importlib.import_module(module_path)
        except ModuleNotFoundError as load_error:
            print('ERROR\nERROR: Could not load "%s". %s\nERROR' % (
                module_path,
                str(load_error)
            ))
            return

        module_classes = inspect.getmembers(module, inspect.isclass)

        for _, plugin_class in module_classes:
            if issubclass(plugin_class, SuperGSLPlugin):
                # Sometimes we encounter a plugin more than once. Ignore already
                # initialized plugins.
                if plugin_class in self._resolved_plugins:
                    continue
                self._resolved_plugins.add(plugin_class)

                print('Registering plugin...', plugin_class)

                plugin_class(self._symbol_table, self._compiler_settings)
