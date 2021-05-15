"""Support for SuperGSL's plugin infrastructure."""

import inspect
import importlib
from typing import Dict
from supergsl.core.exception import ConfigurationError, NotFoundError
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.symbol_table import SymbolTable


class SuperGSLPlugin(object):
    """Base class for defining a SuperGSL Plugin.

    An example of writing a function for SuperGSL:

    class ExamplePlugin(SuperGSLPlugin):

        # Types of plugins
        #    - SuperGSL Function


        def register(self, compiler_settings):
            cut_func = FunctionRegistration(
                'cut', self.cut, return_value=sgsl_types.PART_LIST, help=self.cut.__docstr__)

            cut_fun.add_argument('target_gene', sgsl_types.PART)
            cut_fun.add_argument('target_genome', sgsl_types.GENOME)
            cut_fun.add_argument('target_genome', sgsl_types.GENOME)
            cut_fun.add_argument('num_results', int)

        def cut(self, sgsl_args):
            ""
            Generate gRNA for CRISPR/Cas9

            Example sGSL syntax:
                from example import cut
                cut(HO, CAS9, S288C, results=5)
            ""
    """

    def __init__(self, symbol_table : SymbolTable, compiler_settings : dict):
        self.symbol_table = symbol_table
        self.functions : Dict[str, SuperGSLFunctionDeclaration] = {}
        self.register(compiler_settings)

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Import a identifier and register it in the symbol table."""
        if identifier not in self.functions:
            raise NotFoundError('%s not found in module.' % identifier)

        symbol_table.insert(alias or identifier, self.functions[identifier])


    def register_function(
        self,
        import_path : str,
        function_name : str,
        function_declaration : SuperGSLFunctionDeclaration
    ):
        """Register a function making it available for import in SuperGSL."""
        nested_symbol_table = self.symbol_table.nested_scope('imports')
        self.functions[function_name] = function_declaration
        nested_symbol_table.insert(import_path, self)


    def register(self, compiler_settings : dict):
        """Register Functions, enums, etc that the plugin provides.

        Example: symbol_table.register(import_path, mod_class)
        """
        pass


class PluginProvider(object):
    """Resolve and register SuperGSL plugins."""

    def __init__(self, symbol_table: SymbolTable, compiler_settings : dict):
        self._plugins: Dict[str, SuperGSLPlugin] = {}
        self._symbol_table: SymbolTable = symbol_table
        self._compiler_settings = compiler_settings

        if 'plugins' not in compiler_settings:
            raise ConfigurationError(
                'No plugins have been defined. Check your supergGSL settings.')

        for plugin_path in compiler_settings['plugins']:
            print('Resolving plugin "%s"' % plugin_path)

            self.resolve_plugins_from_config(plugin_path)


    def resolve_plugins_from_config(self, module_path: str) -> None:
        """Attempt to resolve and register a plugin at a specific path."""
        module = importlib.import_module(module_path)
        module_classes = inspect.getmembers(module, inspect.isclass)

        for _, plugin_class in module_classes:
            if issubclass(plugin_class, SuperGSLPlugin):
                print('Registering plugin...', plugin_class)

                plugin_class(self._symbol_table, self._compiler_settings)
