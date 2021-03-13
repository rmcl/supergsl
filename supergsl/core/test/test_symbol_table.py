"""Unit tests for the symbol table."""
import unittest
from mock import Mock
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.provider import SuperGSLProvider

class SymbolTableTestCase(unittest.TestCase):

    def setUp(self):
        self.table = SymbolTable()

    def test_register_and_get_providers(self):
        """Register providers for a specific path and then retrieve it."""

        class MockProvider(SuperGSLProvider):
            pass

        provider1 = MockProvider()
        provider2 = MockProvider()

        self.table.register('impor.this.path', provider1)
        self.table.register('impor.this.path', provider2)

        providers = self.table.get_providers_for_path('impor.this.path')

        self.assertEquals(providers, [
            provider1,
            provider2
        ])

