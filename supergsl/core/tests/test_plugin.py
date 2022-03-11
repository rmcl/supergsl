"""Unit tests for the Compiler Pipeline."""
from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.plugin import PluginProvider, SuperGSLPlugin


class PluginProviderTestCase(TestCase):
    """Test the plugin provider."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        symbol_table = self.fixtures.mk_global_symbol_table()
        self.plugin_provider = PluginProvider(symbol_table, {
            'plugins': []
        })

    def test_get_adhoc_plugin(self):
        """Test that plugin provider instantiates a base plugin and doesn't reinstantiate it."""

        plugin = self.plugin_provider.get_adhoc_plugin()
        self.assertTrue(isinstance(plugin, SuperGSLPlugin))
        plugin2 = self.plugin_provider.get_adhoc_plugin()
        self.assertEqual(plugin, plugin2)
