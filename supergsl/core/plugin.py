import inspect
import importlib
from supergsl.core.config import settings

class SuperGSLPlugin(object):

    def register(self, symbol_table):
        """Register Functions, enums, etc that the plugin provides.

        Example: symbol_table.register(import_path, mod_class)
        """
        pass


class PluginProvider(object):

    def __init__(self, symbol_table):
        self._plugins : Dict[str, SuperGSLPluginConfig] = {}
        self._symbol_table = symbol_table

        if 'plugins' not in settings:
            raise ConfigurationException('No plugins have been defined. Check your supergGSL settings.')

        for plugin_path in settings['plugins']:
            print('Resolving plugin "%s"' % plugin_path)

            self.resolve_plugins_from_config(plugin_path)


    def resolve_plugins_from_config(self, module_path):
        module = importlib.import_module(module_path)
        module_classes = inspect.getmembers(module, inspect.isclass)

        for _, plugin_class in module_classes:
            if issubclass(plugin_class, SuperGSLPlugin):
                print('Registering plugin...', plugin_class)

                plugin_inst = plugin_class()
                plugin_inst.register(self._symbol_table)
