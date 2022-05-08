"""Evaluate a SuperGSL Program."""
from typing import Any, Dict, Optional, Callable, Union, List, cast

from supergsl.core.types import SuperGSLType
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.types.builtin import (
    NucleotideSequence,
    AminoAcidSequence,
    Collection,
    SliceAndInvertCollection,
    SliceInvertMixin
)
from supergsl.utils.resolve import resolve_import
from supergsl.core.sequence import SequenceStore
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    AssemblyLevelDeclaration
)
from supergsl.core.types.position import (
    Position,
    Slice
)

from supergsl.core.constants import (
    UNAMBIGUOUS_DNA_SEQUENCE,
    UNAMBIGUOUS_PROTEIN_SEQUENCE,
    NUMBER_CONSTANT,
    STRING_CONSTANT,
    FIVE_PRIME,
    THREE_PRIME
)

from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.exception import (
    FunctionInvokeError,
    SuperGSLTypeError,
    SymbolNotFoundError,
    FunctionNotFoundError,
    SuperGSLError
)

from supergsl.lang.backend import BackendPipelinePass
from supergsl.lang.ast import (
    Node,
    Program,
    Import,
    Assembly,
    SymbolReference,
    VariableDeclaration,
    ListDeclaration,
    FunctionInvocation,
    Slice as AstSlice,
    SlicePosition as AstSlicePosition,
    SequenceConstant,
    Constant
)
from supergsl.lang.exception import BackendError

#pylint: disable=E1136


class EvaluatePass(BackendPipelinePass):
    """Traverse the AST to execute the GSL Program."""

    def __init__(self, symbol_table : SymbolTable):
        self.symbol_table = symbol_table
        self.sequence_store = cast(
            SequenceStore,
            self.symbol_table.lookup('sequences'))

    def get_node_handlers(self) -> Dict[Optional[str], Callable]:
        """Define method handlers for each node in the AST."""
        return {
            'Program': self.visit_program,
            'Import': self.visit_import,
            'VariableDeclaration': self.visit_variable_declaration,
            'ListDeclaration': self.visit_list_declaration,
            'FunctionInvocation': self.visit_function_invocation,
            'Assembly': self.visit_assembly,
            'SymbolReference': self.visit_symbol_reference,
            'Slice': self.visit_slice,
            'SlicePosition': self.visit_slice_position,
            'Constant': self.visit_constant,
            'SequenceConstant': self.visit_sequence_constant,
        }

    def perform(self, ast_node : Node):
        """Initiate a traversal of the AST."""
        self.visit(ast_node)
        return ast_node

    def visit(self, node, *args, **kwargs):
        """Perform dyanmic dispatch to determine the handler for a node."""
        handlers : Dict[Optional[str], Callable] = self.get_node_handlers()
        node_type : str = type(node).__name__
        handler_method = handlers.get(node_type, None)
        if not handler_method:
            raise BackendError('Handler for node %s not specified.' % node_type)

        return handler_method(node, *args, **kwargs)

    def visit_program(self, program_node : Program):
        for import_node in program_node.imports:
            self.visit(import_node)

        if program_node.definitions:
            for definition_node in program_node.definitions.definitions:
                self.visit(definition_node)


    def visit_import(self, import_node : Import):
        for program_import in import_node.imports:
            resolve_import(
                self.symbol_table,
                import_node.module_path,
                program_import.identifier,
                program_import.alias)


    def visit_assembly(self, assembly : Assembly) -> SuperGSLType:
        """Evaluate the Assembly node by traversing eval'ing all the child symbol references."""

        level_declarations : List[AssemblyLevelDeclaration] = []
        for symbol_reference in assembly.symbol_references:
            part = self.visit(symbol_reference)
            # somehow we need to get this label to AssemblyFactor or something like it
            level_declaration = AssemblyLevelDeclaration(
                part,
                symbol_reference.label)

            # TODO: We need to do type checking here.
            # Ultimately I think these "parts" can be part collections, parts,
            # and nucleotide constants
            #
            #for part in parts:
            #    if not isinstance(part, Part):
            #        raise Exception(
            #            'Type error. Assembly declaration expected a set of parts. '
            #            'Got a "%s"' % part)

            level_declarations.append(level_declaration)

        return AssemblyDeclaration(assembly.label, level_declarations)

    def visit_variable_declaration(self, variable_declaration : VariableDeclaration) -> None:
        """Evaluate by visiting child expression and assigning the result to the symobl table."""
        expression = variable_declaration.value
        expression_result = self.visit(expression)

        if variable_declaration.type_declaration:
            raise NotImplementedError(
                'Variable explicit type declarations are not currently supported.')

        self.symbol_table.insert(
            variable_declaration.identifier,
            expression_result)


    def visit_symbol_reference(self, symbol_reference : SymbolReference) -> SuperGSLType:
        """Visit a SymbolReference node and return its evaluated symbol."""
        symbol = self.symbol_table.lookup(symbol_reference.identifier)
        symbol = symbol.eval()

        if symbol_reference.slice or symbol_reference.invert:
            part_slice = None
            if symbol_reference.slice:
                part_slice = self.visit(symbol_reference.slice)

            if isinstance(symbol, Collection):
                symbol = SliceAndInvertCollection(
                    symbol, part_slice, symbol_reference.invert)

            elif issubclass(type(symbol), SliceInvertMixin):
                if part_slice:
                    symbol = symbol.slice(part_slice)

                if symbol_reference.invert:
                    symbol = symbol.invert()

            else:
                raise Exception(f'{type(symbol)} is not slice or invertable.')

        return symbol


    def visit_slice(self, slice_node : AstSlice):
        start_position = self.visit(slice_node.start)
        end_position = self.visit(slice_node.end)

        return Slice(start_position, end_position)

    def visit_slice_position(self, slice_position : AstSlicePosition):
        """Convert a SlicePosition node into a Position for the given Part."""

        if slice_position.postfix == 'S':
            rel_to = FIVE_PRIME
        elif slice_position.postfix == 'E':
            rel_to = THREE_PRIME
        else:
            raise SuperGSLError('Unknown postfix position. "%s"' % slice_position.postfix)

        return Position(
            index=slice_position.index,
            relative_to=rel_to,
            approximate=slice_position.approximate)


    def visit_list_declaration(self, list_declaration : ListDeclaration) -> Collection:
        """Instantiate a Collection based on details of list declaration."""
        return Collection([
            self.visit(item_node)
            for item_node in list_declaration.item_nodes
        ])

    def visit_constant(self, constant_node : Constant):
        """Create the appropriate constant object based on `Constant` type."""
        if constant_node.constant_type == NUMBER_CONSTANT:
            return int(constant_node.value)

        if constant_node.constant_type == STRING_CONSTANT:
            return constant_node.value

        raise SuperGSLTypeError('Unknown constant type.')


    def visit_sequence_constant(
        self,
        sequence_constant : SequenceConstant
    ) -> Union[NucleotideSequence, AminoAcidSequence]:
        """Return a Sequence Type based on the constant defined SequenceConstant Node."""

        sequence_type = sequence_constant.sequence_type
        sequence_entry = self.sequence_store.add_from_reference(sequence_constant.sequence)
        if sequence_type == UNAMBIGUOUS_PROTEIN_SEQUENCE:
            return AminoAcidSequence(sequence_entry)

        if sequence_type == UNAMBIGUOUS_DNA_SEQUENCE:
            return NucleotideSequence(sequence_entry)

        raise SuperGSLTypeError('Unhandled sequence type "%s"' % sequence_type)

    def visit_function_invocation(self, function_invoke_node : FunctionInvocation) -> SuperGSLType:
        """Evaluate this node by initializing and executing a SuperGSLFunction."""

        try:
            function_declaration = self.symbol_table.lookup(function_invoke_node.identifier)
        except SymbolNotFoundError as error:
            raise FunctionNotFoundError(
                'Function %s has not been imported.' % (
                    function_invoke_node.identifier)) from error

        if not isinstance(function_declaration, SuperGSLFunctionDeclaration):
            raise SuperGSLTypeError(
                '{} is not of type FunctionDeclaration. It is a "{}"'.format(
                    function_invoke_node.identifier,
                    type(function_declaration)))

        function_inst = function_declaration.eval()

        child_arguments = []
        children = function_invoke_node.children
        if children:
            child_arguments = [
                self.visit(child)
                for child in children.definitions
            ]

        params = function_invoke_node.params
        positional_arguments = []
        if params:
            for param_value in params:
                positional_arguments.append(self.visit(param_value))

        return function_inst.evaluate_arguments_and_execute(
            positional_arguments,
            child_arguments)
