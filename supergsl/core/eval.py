"""Define the mechanism of SuperGSLFunction and an AST pass to invoke those functions."""
from typing import List, Dict, Set

from supergsl.core.types import SuperGSLType
from supergsl.core.backend import DepthFirstNodeFilteredPass
from supergsl.core.types.builtin import Collection
from supergsl.utils import display_symbol_table

#pylint: disable=E1136


class DependencyGraph:
    """A graph representing evaluation order dependencies.

    The graph is represented as an adjacency list.

    """

    def __init__(self):
        self.node_edges : Dict[SuperGSLType, Set[SuperGSLType]]= {}

    def add_edge(self, node : SuperGSLType, dependent_node : SuperGSLType) -> None:
        if node not in self.node_edges:
            self.node_edges[node] = set()

        if dependent_node not in self.node_edges:
            self.node_edges[dependent_node] = set()

        self.node_edges[node].add(dependent_node)

    def remove_node(self, node : SuperGSLType) -> None:
        """Remove a node from the graph."""
        del self.node_edges[node]
        for _, adjacent_nodes in self.node_edges.items():
            adjacent_nodes.discard(node)

    def get_available_nodes(self) -> List[SuperGSLType]:
        """Return a list of nodes with no dependencies in the graph"""
        available_nodes : List[SuperGSLType] = []
        for node, adjacent_nodes in self.node_edges.items():
            if len(adjacent_nodes) == 0:
                available_nodes.append(node)

        return available_nodes


class EvaluatePass(DepthFirstNodeFilteredPass):
    """Traverse the AST to build a depenency graph. Then execute all nodes in dependency order."""

    def get_node_handlers(self):
        return {
            'VariableDeclaration': self.visit_variable_declaration_node,
            'ListDeclaration': self.visit_list_declaration_node,
        }

    def before_pass(self, ast):
        self.dependency_graph = DependencyGraph()

        display_symbol_table(self.symbol_table)
        return ast

    def visit_list_declaration_node(self, node):
        #node.collection = Collection(type(node.items[0]))
        for list_item_node in node.items:
            self.dependency_graph.add_edge(node, list_item_node)

        return node
    def visit_variable_declaration_node(self, node):

        # Somehow set / or confirm this variable is defined in the symbol table.

        self.dependency_graph.add_edge(node, node.value)
        return node
