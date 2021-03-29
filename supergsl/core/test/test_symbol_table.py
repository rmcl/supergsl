"""Unit tests for the symbol table."""
import unittest
from mock import Mock
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.exception import SymbolNotFoundError

class SymbolTableTestCase(unittest.TestCase):
    """Testcases to evaluate the SymbolTable class."""

    def test_create_insert_lookup_symbol_table(self):
        """Register providers for a specific path and then retrieve it."""

        mock_symbol_1 = Mock()
        table = SymbolTable('STName', None)
        table.insert('HELLO', mock_symbol_1)

        result = table.lookup('HELLO')
        self.assertEqual(result, mock_symbol_1)

    def test_create_nested_scopes(self):
        """Test that we create a nested scope and it can store same values as parent."""
        mock_symbol_1 = Mock()
        mock_symbol_2 = Mock()

        table = SymbolTable('root', None)
        table.insert('HELLO', mock_symbol_1)

        child_table = table.nested_scope('child_scope')

        self.assertRaises(SymbolNotFoundError, child_table.lookup, 'HELLO')

        child_table.insert('HELLO', mock_symbol_2)

        self.assertEqual(table.lookup('HELLO'), mock_symbol_1)
        self.assertEqual(child_table.lookup('HELLO'), mock_symbol_2)
