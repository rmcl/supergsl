import inspect
import importlib
from supergsl.core.config import settings
from .function import FunctionSymbolTable, SuperGSLFunction


class PluginProvider(object):

    def __init__(self, symbol_registry):
        self._plugins : Dict[str, SuperGSLPluginConfig] = {}

        self.symbol_registry = symbol_registry

        if 'plugins' not in settings:
            raise ConfigurationException('No plugins have been defined. Check your supergGSL settings.')

        for plugin_path in settings['plugins']:
            print('Resolving plugin "%s"' % plugin_path)

            self.resolve_plugins_from_config(plugin_path)


    def resolve_plugins_from_config(self, module_path):
        module = importlib.import_module(module_path)
        module_classes = inspect.getmembers(module, inspect.isclass)

        function_symbol_table = self.symbol_registry.get_table('functions')

        function_defined = False
        for name, mod_class in module_classes:

            mod_class_name = getattr(mod_class, 'name', None)
            if not issubclass(mod_class, SuperGSLFunction) or not mod_class_name:
                continue

            print('Registering plugin...', mod_class)
            function_symbol_table.register_function(mod_class)
            function_defined = True

        #if not function_defined:
        #    raise ConfigurationException('Plugin "%s" did not define anything.' % module_path)
