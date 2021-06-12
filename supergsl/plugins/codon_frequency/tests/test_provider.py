import unittest
from unittest.mock import Mock
from supergsl.core.types.codon import CodonFrequencyTable
from supergsl.plugins.codon_frequency.provider import CodonFrequencyTableProvider

class CodonFrequencyTableProviderTestCase(unittest.TestCase):
    """Test case for CodonFrequencyTableProvider."""

    def setUp(self):
        self.provider = CodonFrequencyTableProvider({})

    def test_get_table(self):
        """Test that get table returns a dictionary of table info."""
        table = self.provider.get_table('s_cerevisiae_4932')
        self.assertTrue(isinstance(table, CodonFrequencyTable))
        self.assertEqual(table.name, 's_cerevisiae_4932')

    def test_resolve_import_inserts_into_symbol_table(self):
        symbol_table = Mock()
        self.provider.get_table = Mock(return_value='HELLO')
        self.provider.resolve_import(
            symbol_table,
            's_cerevisiae_4932',
            'thesacc'
        )

        symbol_table.insert.assert_called_once_with('thesacc', 'HELLO')
