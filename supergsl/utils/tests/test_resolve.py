from unittest import TestCase
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.symbol_table import SymbolTable
from supergsl.plugins.builtin.providers import FastaPartProvider
from supergsl.utils.resolve import resolve_provider_import


class ResolveUtilTestCases(TestCase):
    """Test the helper methods."""

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.symbol_table = self.fixtures.mk_global_symbol_table()

    def test_resolve_provider_import_by_class_str(self):
        """Find provider provided as a str import."""


        provider = resolve_provider_import(
            'supergsl.plugins.builtin.providers.FastaPartProvider',
            self.symbol_table)

        self.assertEqual(provider, FastaPartProvider)

    def test_resolve_provider_import_by_symbol_table(self):
        """Resolve a provider thats already in the symbol table."""

        providers = self.symbol_table.enter_nested_scope('available_imports')
        providers.insert('pasta_fasta', FastaPartProvider)

        provider = resolve_provider_import(
            'pasta_fasta',
            self.symbol_table)

        self.assertEqual(provider, FastaPartProvider)
