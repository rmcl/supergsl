"""Tests for the eval module."""
import unittest
from mock import Mock, call
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.eval import EvaluatePass
from supergsl.core.ast import (
    Program,
    Import,
    ImportIdentifier
)

class EvaluatePassTestCase(unittest.TestCase):
    """Testcases to evaluate the EvaluatePass class."""

    def setUp(self):
        self.symbol_table = SymbolTable('global', None)
        self.import_table = self.symbol_table.nested_scope('imports')
        self.eval_pass = EvaluatePass(self.symbol_table)
        self.eval_pass.visit = Mock()

    def test_visit_program(self):
        """Test visiting a SuperGSL program AST Node."""
        definitions = Mock(
            definitions = ['d1', 'd2'])

        program = Program([1,2,3], definitions)

        self.eval_pass.visit_program(program)

        self.eval_pass.visit.assert_has_calls([
            call(1),
            call(2),
            call(3),
            call('d1'),
            call('d2')
        ])

    def test_visit_import(self):
        """Test visiting a Import AST node."""
        provider_mock = Mock()
        self.import_table.insert('mod_path.here', provider_mock)

        import_node = Import(
            [
                'mod_path',
                'here'
            ], [
                ImportIdentifier('hello', None),
                ImportIdentifier('boop', None)
            ])

        self.eval_pass.visit_import(import_node)

        provider_mock.resolve_import.assert_has_calls([
            call(self.symbol_table, 'hello', None),
            call(self.symbol_table, 'boop', None)
        ])
