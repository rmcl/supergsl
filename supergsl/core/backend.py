class BackendException(Exception):
    pass


class BackendPipelinePass(object):
    name = None

    def get_pass_name(self):
        if not self.name:
            return type(self).__name__

        return self.name

    def perform(self, ast):
        raise NotImplemented('Must subclass and implement perform')


class BreadthFirstNodeFilteredPass(BackendPipelinePass):
    """Perform a breadth first traversal of the AST and only visit a subset of node types."""

    name = None

    def visit(self, node):
        """Visit a node.

        This method is called for every node in the AST. The default implementation uses
        the `get_node_handlers` method to delegate calls to handlers for specific node
        type.

        If you want a catch-all handler specify "None" in get_node_handlers. If you
        just want a simple pass you can override this method.
        """
        handlers = self.get_node_handlers()
        handler_method = None
        node_type = type(node).__name__
        try:
            handler_method = handlers[node_type]
        except KeyError:
            pass

        # Node specific handler defined for this node_type
        # see if there is a default handler defined.
        try:
            handler_method = handlers[None]
        except KeyError:
            pass

        if handler_method:
            handler_method(node)

    def get_node_handlers(self):
        """Define handler methods for each node type.

        Returns a dictionary of
        """
        return {}

    def perform(self, ast):
        node_visit_queue = ast.child_nodes().copy()

        while len(node_visit_queue) > 0:
            cur_node = node_visit_queue.pop(0)

            self.visit(cur_node)
            node_visit_queue += cur_node.child_nodes()

        return ast
