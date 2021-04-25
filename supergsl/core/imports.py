from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.parts import Part, LazyLoadedPart
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.exception import SuperGSLTypeError
from supergsl.core.types import resolve_type


class ResolveImportsPass(BreadthFirstNodeFilteredPass):
    """Traverse the AST and attempt to resolve imported symbols."""

    def get_node_handlers(self):
        return {
            'Import': self.visit_import_node,
            'SymbolReference': self.visit_symbol_reference_node,
            'FunctionInvocation': self.visit_function_invoke_node,
            'VariableDeclaration': self.visit_variable_declaration_node,
        }

    def before_pass(self, ast):
        return ast

    def visit_import_node(self, node):
        """Visit AST import nodes and resolve them via the registered providers."""
        import_table = self.symbol_table.nested_scope('imports')

        for program_import in node.imports:
            module_path = '.'.join(node.module)

            provider = import_table.lookup(module_path)
            node.symbol = provider.resolve_import(
                self.symbol_table,
                program_import.identifier,
                program_import.alias)

        return node

    def visit_variable_declaration_node(self, node):
        """

        node.variable = self.symbol_table.resolve_symbol(
            VARIABLE_MODULE_PATH,
            node.identifier)

        supergsl_type = resolve_type(node.type_declaration_node.identifier)

        if node.type_declaration:
            node.variable.set_type(supergsl_type)

        if node.type_declaration:
            node.variable.set_value(node.value)


        """
        return node

    def visit_function_invoke_node(self, node):
        function_declaration = self.symbol_table.lookup(node.identifier)
        if not isinstance(function_declaration, SuperGSLFunctionDeclaration):
            raise SuperGSLTypeError('{} is not of type Function. It is a "{}"'.format(
                node.identifier,
                type(function_declaration)
            ))

        node.set_function_declaration(function_declaration)

        ## TODO: DO TYPE CHECKING TO ASSERT THIS SYMBOL IS A FUNCTION
        #node.expected_type = node.function.get_return_type()
        return node

    def visit_symbol_reference_node(self, node):
        part_symbol = self.symbol_table.lookup(node.identifier)
        if isinstance(part_symbol, LazyLoadedPart):
            part_symbol = part_symbol.instantiate()

        if not isinstance(part_symbol, Part):
            raise SuperGSLTypeError('{} is not of type Part. It is a "{}"'.format(
                node.identifier,
                type(part_symbol)
            ))

        node.set_referenced_object(part_symbol)
        return node
