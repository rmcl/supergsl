"""Tests for the eval module."""
import unittest
from mock import Mock, call, patch
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.eval import EvaluatePass
from supergsl.core.constants import (
    STRING_CONSTANT,
    NUMBER_CONSTANT
)
from supergsl.core.ast import (
    Program,
    Import,
    ImportIdentifier,
    Assembly,
    VariableDeclaration,
    SymbolReference,
    Slice,
    Constant
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

    def test_visit_assembly(self):
        """Visit assembly visits each part and builds an AssemblyDeclaration."""
        self.eval_pass.visit.return_value = 'BOOM'
        parts = [
            Mock(),
            Mock()
        ]
        assembly_node = Assembly(parts, 'LABEL1112')

        assembly_declaration = self.eval_pass.visit_assembly(assembly_node)
        self.assertEqual(assembly_declaration.get_label(), 'LABEL1112')
        self.assertEqual(assembly_declaration.get_parts(), [
            'BOOM',
            'BOOM'
        ])
        self.eval_pass.visit.assert_has_calls([
            call(parts[0]),
            call(parts[1])
        ])

    def test_visit_variable_declaration(self):
        """Variable declaration visits the expression and inserts result into symbol table."""
        self.eval_pass.visit.return_value = 'RESULT'
        var_declare = VariableDeclaration('IDENT', None, 'BOOMVAL')

        self.eval_pass.visit_variable_declaration(var_declare)

        self.eval_pass.visit.assert_has_calls([
            call('BOOMVAL')
        ])

        actual_result = self.symbol_table.lookup('IDENT')
        self.assertEqual('RESULT', actual_result)

    def test_visit_symbol_reference(self):
        """Symbol Reference visit should retrieve a symbol from the symbol table."""
        supergsl_type = Mock()
        supergsl_type.eval.return_value = 'YES!'

        self.symbol_table.insert('IDENT', supergsl_type)
        symbol_ref = SymbolReference('IDENT', None, False)

        result = self.eval_pass.visit_symbol_reference(symbol_ref)

        self.assertEqual(result, 'YES!')

    def test_visit_slice(self):
        """Visit Slice should visit start and end positions and create a child part."""
        start = Mock()
        end = Mock()
        parent_part = Mock(identifier='HIII')
        expected_child_part = Mock()
        parent_part.get_child_part_by_slice.return_value = expected_child_part
        self.eval_pass.visit.return_value = 'VISIT-RETURN-VAL'

        slice_node = Slice(start, end)
        slice_node.get_slice_str = Mock(return_value='poop')

        new_part = self.eval_pass.visit_slice(slice_node, parent_part)

        self.assertEqual(new_part, expected_child_part)
        self.eval_pass.visit.assert_has_calls([
            call(start, parent_part),
            call(end, parent_part)
        ])

        parent_part.get_child_part_by_slice.assert_called_once_with(
            'HIII[poop]',
            'VISIT-RETURN-VAL',
            'VISIT-RETURN-VAL'
        )

    def test_visit_slice_position(self):
        """Visit slice position should convert a slice position to a SeqPosition"""

        parent_part = Mock()
        slice_position = Mock()

        convert_util_path = 'supergsl.core.eval.convert_slice_position_to_seq_position'
        with patch(convert_util_path) as convert_patch:
            convert_patch.return_value = 'HELLOOO'

            result = self.eval_pass.visit_slice_position(
                slice_position, parent_part)

            convert_patch.assert_called_once_with(parent_part, slice_position)
            self.assertEqual(result, 'HELLOOO')

    # TODO: visit_list_declaration

    def test_visit_constant_number(self):
        """Visit a constant node resolves to a number."""
        constant_node = Constant(55, NUMBER_CONSTANT)

        result = self.eval_pass.visit_constant(constant_node)
        self.assertEqual(result, 55)

    def test_visit_constant_string(self):
        """Visit a constant node resolves to a string."""
        constant_node = Constant('party', STRING_CONSTANT)

        result = self.eval_pass.visit_constant(constant_node)
        self.assertEqual(result, 'party')

    # TODO: visit_sequence_constant
