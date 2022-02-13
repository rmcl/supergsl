from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.symbol_table import SymbolTable
from supergsl.plugins.builtin.providers import FastaPartProvider
from supergsl.utils.resolve import resolve_provider_import


class ResolveUtilTestCases(TestCase):
    """Test the helper methods."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()

    def test_find_or_import_provider(self):
        """Find provider provided as a str import."""


        provider = resolve_provider_import(
            'supergsl.plugins.builtin.providers.FastaPartProvider',
            self.fixtures.mk_global_symbol_table())

        self.assertEqual(provider, FastaPartProvider)
