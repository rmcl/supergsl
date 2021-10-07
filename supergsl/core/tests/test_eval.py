"""Tests for the eval module."""
import unittest
from unittest.mock import Mock, call, patch
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.eval import EvaluatePass
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.constants import (
    THREE_PRIME,
    STRING_CONSTANT,
    NUMBER_CONSTANT,
    UNAMBIGUOUS_DNA_SEQUENCE,
    UNAMBIGUOUS_PROTEIN_SEQUENCE
)
from supergsl.core.ast import (
    Program,
    Import,
    ImportIdentifier,
    Assembly,
    VariableDeclaration,
    SymbolReference,
    Slice,
    SlicePosition,
    Constant,
    ListDeclaration,
    SequenceConstant,
    FunctionInvocation
)

from supergsl.core.types.slice import Position
from supergsl.core.types.builtin import (
    Collection,
    NucleotideSequence,
    AminoAcidSequence
)

from supergsl.core.exception import (
    SuperGSLTypeError,
    FunctionNotFoundError
)

class EvaluatePassTestCase(unittest.TestCase):
    """Testcases to evaluate the EvaluatePass class."""

    def setUp(self):
        self.symbol_table = SymbolTable('global', None)
        self.import_table = self.symbol_table.enter_nested_scope('imports')
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
        provider_mock = Mock(__class__ = SuperGSLProvider)
        self.import_table.insert('mod_path.here', provider_mock)
        provider_mock.resolve_import.return_value = {}

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
            call('hello', None),
            call('boop', None)
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
        self.assertEqual(assembly_declaration.label, 'LABEL1112')
        self.assertEqual(assembly_declaration.get_levels_by_factor_type('Part'), {
            'BOOM',
        })
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
        """Visit Slice should visit start and end positions and create a Slice type."""
        start = Mock()
        end = Mock()

        self.eval_pass.visit.return_value = position_mock = Position(0)

        slice_node = Slice(start, end)
        slice_obj = self.eval_pass.visit_slice(slice_node)

        self.assertEqual(slice_obj.start, position_mock)
        self.assertEqual(slice_obj.end, position_mock)

        self.eval_pass.visit.assert_has_calls([
            call(start),
            call(end)
        ])

    def test_visit_slice_position(self):
        """Visit slice position should convert a SlicePosition AST node to a Position"""

        slice_position = SlicePosition(123, 'E', False)

        result = self.eval_pass.visit_slice_position(
            slice_position)

        self.assertEqual(result.index, 123)
        self.assertEqual(result.relative_to, THREE_PRIME)
        self.assertEqual(result.approximate, False)

    def test_visit_list_declaration(self):
        """Create a collection with result of visiting each item node in declaration."""
        item_nodes = [
            Mock(),
            Mock(),
            Mock(),
        ]
        list_declare_node = ListDeclaration(item_nodes)

        result = self.eval_pass.visit_list_declaration(list_declare_node)

        self.assertEqual(type(result), Collection)
        self.eval_pass.visit.assert_has_calls([
            call(item_node)
            for item_node in item_nodes
        ])

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

    def test_visit_dna_sequence_constant(self):
        """Create a GSL Sequence Type object based on the declared sequence constant."""
        constant_dna_node = SequenceConstant(
            'ATGC', UNAMBIGUOUS_DNA_SEQUENCE)

        result = self.eval_pass.visit_sequence_constant(constant_dna_node)
        self.assertEqual(type(result), NucleotideSequence)
        self.assertEqual(result.sequence, 'ATGC')

    def test_visit_protein_sequence_constant(self):
        """Create the right type for a protein sequence constant."""
        constant_protein_node = SequenceConstant(
            'MATTTGAC*', UNAMBIGUOUS_PROTEIN_SEQUENCE)

        result = self.eval_pass.visit_sequence_constant(constant_protein_node)
        self.assertEqual(type(result), AminoAcidSequence)
        self.assertEqual(result.sequence, 'MATTTGAC*')

    def test_visit_sequence_constant_weird_type(self):
        """Raise exception for a unknown type of constant."""
        constant_protein_node = SequenceConstant(
            'MATTTGAC*', 'PARTY_TYPE')

        self.assertRaises(
            SuperGSLTypeError,
            self.eval_pass.visit_sequence_constant,
            constant_protein_node)

    def test_visit_function_invocation_function_not_defined(self):
        """Raise exception if the desired function is not in the symbol table."""
        function_invoke_node = FunctionInvocation(
            'great_function',
            Mock(),
            [],
            None)

        self.assertRaises(
            FunctionNotFoundError,
            self.eval_pass.visit_function_invocation,
            function_invoke_node)

    def test_visit_function_invocation_not_function_declaration(self):
        """Raise exception if symbol in symbol table is not a function declare."""
        self.symbol_table.insert('great_function', AminoAcidSequence('MAA*'))
        function_invoke_node = FunctionInvocation(
            'great_function',
            Mock(),
            [],
            None)

        self.assertRaises(
            SuperGSLTypeError,
            self.eval_pass.visit_function_invocation,
            function_invoke_node)

    def test_visit_function_invocation(self):
        function_declare = SuperGSLFunctionDeclaration(Mock(), {})
        func_inst = Mock()
        function_declare.eval = Mock(return_value = func_inst)
        func_inst.evaluate_arguments_and_execute.return_value = 'HELLO world!'
        self.symbol_table.insert('great_function', function_declare)

        function_invoke_node = FunctionInvocation(
            'great_function', [], [], None)

        result = self.eval_pass.visit_function_invocation(function_invoke_node)
        self.assertEqual(result, 'HELLO world!')
        func_inst.evaluate_arguments_and_execute.assert_called_once_with(
            [], [])
