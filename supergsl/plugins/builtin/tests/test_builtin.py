"""Unit tests for the Builtin Plugin."""
import unittest
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin import BuiltinPlugin

class BuiltinPluginTestCase(unittest.TestCase):
    """Test for the builtin plugin."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_initialize_builtin_plugin_functions(self):
        """Confirm the plugin registers expected functions."""

        symbol_table = self.fixtures.mk_global_symbol_table()
        builtin_plugin = BuiltinPlugin(symbol_table, {})

        self.assertEqual(list(builtin_plugin.functions.keys()), [
            'fuse',
            'synthetic_oligos',
            'detail',
            'help',
            'declare',
            'output_json',
            'output_sbol',
            'output_genbank'
        ])
