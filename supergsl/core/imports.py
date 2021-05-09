from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.exception import SuperGSLTypeError


class ResolveImportsPass(BreadthFirstNodeFilteredPass):
    """Traverse the AST and attempt to resolve imported symbols."""

    def get_node_handlers(self):
        return {
            'Import': self.visit_import_node,
            'VariableDeclaration': self.visit_variable_declaration_node,
            'SymbolReference': self.visit_symbol_reference_node,
            'FunctionInvocation': self.visit_function_invoke_node,
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


    def visit_function_invoke_node(self, node):
        """Visit AST function invoke node and lookup the corresponding function."""
        function_declaration = self.symbol_table.lookup(node.identifier)
        if not isinstance(function_declaration, SuperGSLFunctionDeclaration):
            raise SuperGSLTypeError('{} is not of type Function. It is a "{}"'.format(
                node.identifier,
                type(function_declaration)
            ))

        node.set_function_declaration(function_declaration)

        # TODO: DO TYPE CHECKING TO ASSERT THIS SYMBOL IS A FUNCTION
        # node.expected_type = node.function.get_return_type()

        return node

    def visit_variable_declaration_node(self, node):
        """Visit each variable declaration and define the variable."""
        node.set_table_reference(self.symbol_table, node.identifier)

        # Set a value to the symbol table so we know the variable is declared.
        self.symbol_table.insert(node.identifier, None)
        return node

    def visit_symbol_reference_node(self, node):
        """Visit an AST SymbolReference node and lookup the corresponding symbol."""

        # Lookup the symbol to see if it has been defined.
        self.symbol_table.lookup(node.identifier)

        node.set_table_reference(self.symbol_table, node.identifier)
        return node
