"""Implement an AST Traversal used by tests to introspect the output of the compiler"""

from supergsl.core.backend import BreadthFirstNodeFilteredPass


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
        self.parts = []
        self.assemblies = []

    def visit_symbol_reference_node(self, node):
        # Todo(rmcl): make this more generic to handle symbol references other
        # than parts.
        self.parts.append(node.part)

    def visit_assembly_node(self, node):
        assembly_parts = [
            part.identifier
            for part in node.parts
        ]

        assembly_idx = len(self.assemblies)
        self.assemblies.append({
            'identifier': assembly_idx,
            'parts': assembly_parts
        })

    def get_assemblies(self):
        return self.assemblies

    def get_parts(self):
        return self.parts
