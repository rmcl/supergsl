"""Support for SuperGSL's plugin infrastructure."""

import inspect
import importlib
from typing import Dict
from supergsl.core.exception import ConfigurationError
from supergsl.core.symbol_table import SymbolTable


class SuperGSLPlugin(object):
    """Base class for defining a SuperGSL Plugin.

    An example of writing a function for SuperGSL:

    class ExamplePlugin(SuperGSLPlugin):

        # Types of plugins
        #    - SuperGSL Function


        def register(self, context):
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

    def register(self, symbol_table, compiler_settings):
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

                plugin_inst = plugin_class()
                plugin_inst.register(self._symbol_table, self._compiler_settings)
