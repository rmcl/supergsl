from typing import List, Tuple, Optional, Dict, Callable
from .ast import Node

from supergsl.core.symbol_table import SymbolTable
from supergsl.core.exception import BackendError

#pylint: disable=E1136

# Define a mypy type alias for node handler methods.
ASTNodeHandlerMethod = Callable[[Node], Node]


class BackendPipelinePass(object):
    """Base class for implementing a traversal of the AST."""
    name : Optional[str] = None

    def __init__(self, symbol_table, allow_modification : bool = True):
        self.symbol_table = symbol_table
        self.allow_modification = allow_modification

    def get_pass_name(self):
        if not self.name:
            return type(self).__name__

        return self.name

    def perform(self, ast : Node) -> Node:
        raise NotImplementedError('Must subclass and implement perform')


class BreadthFirstNodeFilteredPass(BackendPipelinePass):
    """Perform a preorder breadth first traversal of the AST and only visit a subset of node types."""

    def visit(self, node : Node, parent_node : Optional[Node]) -> None:
        """Visit a node.

        This method is called for every node in the AST. The default implementation uses
        the `get_node_handlers` method to delegate calls to handlers for specific node
        type.

        If you want a catch-all handler specify "None" in get_node_handlers. If you
        just want a simple pass you can override this method.
        """
        if not node:
            raise BackendError('Past "{}" was passed a null AST.'.format(self.get_pass_name()))

        handler_method = self.get_handler_for_node(node)
        if handler_method:
            result_node = self.call_handler_and_check_result(handler_method, node)
            if result_node != node:
                if not parent_node:
                    raise BackendError(
                        'You cannot update the root Program AST node. Tried to update "%s"' % node)

                parent_node.replace_child_node(node, result_node)

    def get_handler_for_node(self, node : Node) -> Optional[ASTNodeHandlerMethod]:
        """Get the handler for the given AST node.

        If no handler specified, try to return the default node handler. If no
        default handler speciifed then return None.
        """
        handlers = self.get_node_handlers()

        node_type = type(node).__name__
        try:
            return handlers[node_type]
        except KeyError:
            # No Node specific handler defined for this node_type
            # see if there is a default handler defined.
            return handlers.get(None, None)

    def get_node_handlers(self) -> Dict[Optional[str], ASTNodeHandlerMethod]:
        """Define handler methods for each node type.

        Returns a dictionary of
        """
        return {}

    def before_pass(self, ast : Node) -> Node:
        return ast

    def after_pass(self, ast : Node) -> Node:
        return ast

    def call_handler_and_check_result(self, handler, ast_node):
        """Call handler method on an AST node and check to if the node was modified.

        Apply rules based on the initialization of the pipeline.
        """
        new_ast_node = handler(ast_node)

        if not self.allow_modification:
            if new_ast_node:
                raise BackendError(
                    'Handler of "%s" should not modify the AST. Your node handlers should return `None`.' % (
                        handler
                    ))
        else:
            if not new_ast_node:
                raise BackendError(
                    'Handler "%s" pass did not return an AST node object.' % handler)

        return new_ast_node or ast_node

    def perform(self, ast : Node) -> Node:
        ast = self.call_handler_and_check_result(self.before_pass, ast)

        node_visit_queue : List[Tuple[Node, Optional[Node]]] = [(ast, None)]
        while len(node_visit_queue) > 0:
            cur_node, cur_node_parent = node_visit_queue.pop(0)

            self.visit(cur_node, cur_node_parent)

            child_nodes = cur_node.child_nodes()
            if child_nodes:
                for child in cur_node.child_nodes():
                    node_visit_queue.append((child, cur_node))

        ast = self.call_handler_and_check_result(self.after_pass, ast)
        return ast


class DepthFirstNodeFilteredPass(BreadthFirstNodeFilteredPass):
    """Perform a postorder depth first traversal of the AST and only visit a subset of node types."""

    def perform(self, ast : Node) -> Node:
        ast = self.call_handler_and_check_result(self.before_pass, ast)

        # Store both current node and parent.
        node_stack : List[Tuple[Node, Optional[Node]]] = [(ast, None)]
        discovered = set()

        while len(node_stack) > 0:
            cur_node, cur_node_parent = node_stack[-1] # peek

            if cur_node in discovered:
                node_stack.pop()
                self.visit(cur_node, cur_node_parent)
            else:
                discovered.add(cur_node)
                child_nodes = cur_node.child_nodes()
                if child_nodes:
                    node_stack.extend([
                        (child_node, cur_node)
                        for child_node in child_nodes
                    ])

        ast = self.call_handler_and_check_result(self.after_pass, ast)
        return ast
