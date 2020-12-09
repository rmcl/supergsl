from typing import Optional, Dict, Callable
from .ast import Node

from supergsl.core.ast import SymbolRepository


# Define a mypy type alias for node handler methods.
ASTNodeHandlerMethod = Callable[[Node], Node]

class BackendException(Exception):
    pass


class BackendPipelinePass(object):
    """Base class for implementing a traversal of the AST."""
    name : Optional[str] = None

    def __init__(self, symbol_registry):
        self.symbol_registry = symbol_registry

    def get_pass_name(self):
        if not self.name:
            return type(self).__name__

        return self.name

    def perform(self, ast : Node) -> Node:
        raise NotImplementedError('Must subclass and implement perform')


class BreadthFirstNodeFilteredPass(BackendPipelinePass):
    """Perform a breadth first traversal of the AST and only visit a subset of node types."""

    def visit(self, node : Node) -> None:
        """Visit a node.

        This method is called for every node in the AST. The default implementation uses
        the `get_node_handlers` method to delegate calls to handlers for specific node
        type.

        If you want a catch-all handler specify "None" in get_node_handlers. If you
        just want a simple pass you can override this method.
        """
        handlers : Dict[Optional[str], ASTNodeHandlerMethod] = self.get_node_handlers()
        handler_method : Optional[ASTNodeHandlerMethod] = None

        node_type : str = type(node).__name__
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

    def get_node_handlers(self) -> Dict[Optional[str], ASTNodeHandlerMethod]:
        """Define handler methods for each node type.

        Returns a dictionary of
        """
        return {}

    def before_pass(self, ast : Node) -> Node:
        return ast

    def after_pass(self, ast : Node) -> Node:
        return ast

    def perform(self, ast : Node) -> Node:
        ast = self.before_pass(ast)

        if not ast:
            raise BackendException('before_pass of "%s" did not return an AST node object.' % self)

        node_visit_queue = [ast]
        while len(node_visit_queue) > 0:
            cur_node = node_visit_queue.pop(0)

            self.visit(cur_node)
            node_visit_queue += cur_node.child_nodes()

        ast = self.after_pass(ast)
        if not ast:
            raise BackendException('after_pass of "%s" did not return an AST node object.' % self)

        return ast


class DepthFirstNodeFilteredPass(BreadthFirstNodeFilteredPass):
    """Perform a depth first traversal of the AST and only visit a subset of node types."""

    def perform(self, ast : Node) -> Node:
        ast = self.before_pass(ast)

        if not ast:
            raise BackendException('before_pass of "%s" did not return an AST node object.' % self)

        node_stack = [ast]
        discovered = set()

        while len(node_stack) > 0:
            cur_node = node_stack[-1] # peek

            if cur_node in discovered:
                node_stack.pop()
                self.visit(cur_node)
            else:
                discovered.add(cur_node)
                node_stack += cur_node.child_nodes()

        ast = self.after_pass(ast)
        if not ast:
            raise BackendException('after_pass of "%s" did not return an AST node object.' % self)

        return ast
