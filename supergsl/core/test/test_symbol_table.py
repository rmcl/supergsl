import unittest
import mock
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.provider import SuperGSLProvider

class SymbolTableTestCase(unittest.TestCase):

    def setUp(self):
        self.table = SymbolTable()

    def test_register_and_get_providers(self):
        """Register providers for a specific path and then retrieve it."""

        class MockProvider(SuperGSLProvider):
            pass

        class MockProvider2(SuperGSLProvider):
            pass

        self.table.register('impor.this.path', MockProvider)
        self.table.register('impor.this.path', MockProvider2)

        providers = self.table.get_providers_for_path('impor.this.path')

        self.assertEquals(providers, [
            MockProvider,
            MockProvider2
        ])
