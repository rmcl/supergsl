"""Evaluate a SuperGSL Program."""
from typing import Any, Dict, Optional, Callable

from supergsl.core.types import SuperGSLType
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.backend import BackendPipelinePass
from supergsl.core.types.builtin import Collection
from supergsl.core.types.assembly import AssemblyDeclaration
from supergsl.utils import display_symbol_table
from supergsl.core.ast import (
    Node,
    Program,
    Import,
    Assembly,
    SymbolReference,
    VariableDeclaration,
    ListDeclaration,
    FunctionInvocation
)

from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.exception import FunctionInvokeError, SuperGSLTypeError

#pylint: disable=E1136


class EvaluatePass(BackendPipelinePass):
    """Traverse the AST to execute the GSL Program."""

    def __init__(self, symbol_table : SymbolTable):
        self.symbol_table = symbol_table

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
        }

    def perform(self, ast_node : Node):
        """Initiate a traversal of the AST."""
        self.visit(ast_node)

    def visit(self, node):
        """Perform dyanmic dispatch to determine the handler for a node."""
        handlers : Dict[Optional[str], Callable] = self.get_node_handlers()
        node_type : str = type(node).__name__
        handler_method = handlers.get(node_type, None)
        if not handler_method:
            raise Exception('Handler for node %s not specified.' % node_type)

        print('VISIT', node)
        return handler_method(node)

    def visit_program(self, program_node : Program):
        for import_node in program_node.imports:
            self.visit(import_node)

        for definition_node in program_node.definitions.definitions:
            self.visit(definition_node)


    def visit_import(self, import_node : Import):
        import_table = self.symbol_table.nested_scope('imports')

        for program_import in import_node.imports:
            module_path = '.'.join(import_node.module)

            provider = import_table.lookup(module_path)
            provider.resolve_import(
                self.symbol_table,
                program_import.identifier,
                program_import.alias)


    def visit_assembly(self, assembly : Assembly) -> SuperGSLType:
        """Evaluate the Assembly node by traversing eval'ing all the child symbol references."""

        parts = []
        for symbol_reference in assembly.symbol_references:
            parts.append(self.visit(symbol_reference))

        # Todo: We need to do type checking here. Figure out how to move this
        # logic out of AST!!!
        #for part in parts:
        #    if not isinstance(part, Part):
        #        raise Exception(
        #            'Type error. Assembly declaration expected a set of parts. '
        #            'Got a "%s"' % part)

        return AssemblyDeclaration(assembly.label, parts)

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
        symbol = self.symbol_table.lookup(symbol_reference.identifier)
        return symbol.eval(self)

    def visit_list_declaration(self, list_declaration : ListDeclaration) -> Collection:
        return Collection([
            self.visit(item_node)
            for item_node in list_declaration.item_nodes
        ])

    def visit_function_invocation(self, function_invoke_node : FunctionInvocation) -> SuperGSLType:
        """Evaluate this node by initializing and executing a SuperGSLFunction."""

        function_declaration = self.symbol_table.lookup(function_invoke_node.identifier)
        if not isinstance(function_declaration, SuperGSLFunctionDeclaration):
            raise SuperGSLTypeError('{} is not of type Function. It is a "{}"'.format(
                function_invoke_node.identifier,
                type(function_declaration)
            ))

        function_inst = function_declaration.eval(self)

        expected_return_type = function_inst.get_return_type()
        print('expected', expected_return_type)

        eval_params : Dict[Any, Any] = {
            'children': []
        }

        params = function_invoke_node.params
        children = function_invoke_node.children
        if params:
            # TODO: Right now we are using positional indexes for arguments in a
            # rather inelegant way. Reflect on this and attempt to improve.
            for idx in range(len(params)):
                print(params[idx])
                eval_params[idx] = self.visit(params[idx])

        if children:
            eval_params['children'] = [
                self.visit(child)
                for child in children.definitions
            ]

        function_result = function_inst.execute(eval_params)
        if not isinstance(function_result, expected_return_type):
            raise FunctionInvokeError(
                '"%s" Return type does not match expectation. Expected: "%s", Actual: "%s"' % (
                    function_inst,
                    expected_return_type,
                    type(function_result)
                ))

        return function_result
