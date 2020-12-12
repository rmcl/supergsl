from typing import List, Tuple, Optional, Dict, Callable
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
    """Perform a preorder breadth first traversal of the AST and only visit a subset of node types."""

    def visit(self, node : Node, parent_node : Optional[Node]) -> None:
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
        handler_method = handlers.get(node_type, None)
        if not handler_method:
            # No Node specific handler defined for this node_type
            # see if there is a default handler defined.
            handler_method = handlers.get(None, None)

        if handler_method:
            result_node = handler_method(node)

            if result_node is None:
                raise Exception('Node handler "%s" must return either the supplied AST node or a replacement node.' % handler_method)

            print('VISIT', result_node, node, parent_node)
            if result_node != node:
                if not parent_node:
                    raise Exception('You cannot update the root Program AST node. Tried to update "%s"' % node)

                parent_node.replace_child_node(node, result_node)

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

        node_visit_queue : List[Tuple[Node, Optional[Node]]] = [(ast, None)]
        while len(node_visit_queue) > 0:
            cur_node, cur_node_parent = node_visit_queue.pop(0)

            self.visit(cur_node, cur_node_parent)
            node_visit_queue += [
                (child, cur_node)
                for child in cur_node.child_nodes()
            ]

        ast = self.after_pass(ast)
        if not ast:
            raise BackendException('after_pass of "%s" did not return an AST node object.' % self)

        return ast


class DepthFirstNodeFilteredPass(BreadthFirstNodeFilteredPass):
    """Perform a postorder depth first traversal of the AST and only visit a subset of node types."""

    def perform(self, ast : Node) -> Node:
        ast = self.before_pass(ast)

        if not ast:
            raise BackendException('before_pass of "%s" did not return an AST node object.' % self)

        node_stack : List[Tuple[Node, Optional[Node]]] = [(ast, None)]
        discovered = set()

        while len(node_stack) > 0:
            cur_node, cur_node_parent = node_stack[-1] # peek

            if cur_node in discovered:
                node_stack.pop()
                self.visit(cur_node, cur_node_parent)
            else:
                discovered.add(cur_node)
                node_stack += [
                    (child_node, cur_node) # Store both current node and parent.
                    for child_node in cur_node.child_nodes()
                ]

        ast = self.after_pass(ast)
        if not ast:
            raise BackendException('after_pass of "%s" did not return an AST node object.' % self)

        return ast
