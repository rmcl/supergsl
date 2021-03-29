from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.parts import Part, LazyLoadedPart
from supergsl.core.function import SuperGSLFunction
from supergsl.core.exception import SuperGSLTypeError


class ResolveImportsPass(BreadthFirstNodeFilteredPass):
    """Traverse the AST and attempt to resolve imported symbols."""

    def get_node_handlers(self):
        return {
            'Import': self.visit_import_node,
            'Part': self.visit_part_node,
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
        node.function = self.symbol_table.lookup(node.identifier)
        if not isinstance(node.function, SuperGSLFunction):
            raise SuperGSLTypeError('{} is not of type Function. It is a "{}"'.format(
                node.identifier,
                type(node.function)
            ))

        ## TODO: DO TYPE CHECKING TO ASSERT THIS SYMBOL IS A FUNCTION
        #node.expected_type = node.function.get_return_type()
        return node

    def visit_part_node(self, node):

        part_symbol = self.symbol_table.lookup(node.identifier)
        if isinstance(part_symbol, LazyLoadedPart):
            part_symbol = part_symbol.instantiate()

        if not isinstance(part_symbol, Part):
            raise SuperGSLTypeError('{} is not of type Part. It is a "{}"'.format(
                node.identifier,
                type(part_symbol)
            ))

        node.part = part_symbol
        return node
