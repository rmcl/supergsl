"""Implement an AST Traversal used by tests to introspect the output of the compiler"""

from supergsl.lang.backend import BreadthFirstNodeFilteredPass


class TestOutputAstPass(BreadthFirstNodeFilteredPass):
    """AST Traversal used by Integration tests to introspect the output of the compiler."""
    name = 'test'

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly_node,
            'SymbolReference': self.visit_symbol_reference_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        pass

    def visit_symbol_reference_node(self, node):
        pass

    def visit_assembly_node(self, node):
        pass
