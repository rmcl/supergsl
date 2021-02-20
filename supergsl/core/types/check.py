from supergsl.core.backend import DepthFirstNodeFilteredPass


class TypeCheckPass(DepthFirstNodeFilteredPass):
    """Traverse the AST and validate type constraints."""

    def get_node_handlers(self):
        return {
            'FunctionInvocation': self.visit_function_invoke_node,
        }

    #def visit_list_declaration_node(self, node):
    #    """Visit ListDeclaration nodes and confirm all items are of the same type."""
    #    pass


    def visit_function_invoke_node(self, node):
        """Visit FunctionInvocation node and determine the expected return type."""

        node.expected_type = node.function.get_return_type()

        return node
