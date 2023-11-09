"""Unittests for the backend.py module."""
from unittest import TestCase
from unittest.mock import Mock, call
from typing import Optional, List
from supergsl.lang.ast import Node
from supergsl.lang.backend import (
    BreadthFirstNodeFilteredPass,
    DepthFirstNodeFilteredPass
)


class ExampleNode(Node):
    """An example AST node which makes it easy to set child nodes."""
    def __init__(self, name, child_nodes : Optional[List[Node]] = None):
        self.name = name
        self._child_nodes = child_nodes

    def child_nodes(self) -> List[Node]:
        return self._child_nodes or []

    def __repr__(self) -> str:
        return self.name

class OrderedTraversalFilteredPassTestCase(TestCase):
    """Test that depth first traversal works as expected."""

    def setUp(self):
        self.symbol_table = Mock()

        self.leaf1 = ExampleNode('leaf1', [])
        self.leaf2 = ExampleNode('leaf2', [])
        self.mid1 = ExampleNode('mid1', [self.leaf1])
        self.mid2 = ExampleNode('mid2', [self.leaf2])
        self.root = ExampleNode('root', [self.mid1, self.mid2])

    def test_breadth_first_node_traversal_perform(self):
        """Test that the breadth first node traversl perform calls nodes in right order."""

        filter_pass = BreadthFirstNodeFilteredPass(self.symbol_table, False)
        filter_pass.before_pass = Mock(return_value = None)
        filter_pass.after_pass = Mock(return_value = None)
        filter_pass.visit = Mock()

        filter_pass.perform(self.root)

        self.assertEqual(filter_pass.visit.call_args_list, [
            call(self.root, None),
            call(self.mid1, self.root),
            call(self.mid2, self.root),
            call(self.leaf1, self.mid1),
            call(self.leaf2, self.mid2),
        ])

    def test_depth_first_node_traversal_perform(self):
        """Test that the breadth first node traversl perform calls nodes in right order."""

        filter_pass = DepthFirstNodeFilteredPass(self.symbol_table, False)
        filter_pass.before_pass = Mock(return_value = None)
        filter_pass.after_pass = Mock(return_value = None)
        filter_pass.visit = Mock()

        filter_pass.perform(self.root)

        self.assertEqual(filter_pass.visit.call_args_list, [
            call(self.leaf2, self.mid2),
            call(self.mid2, self.root),
            call(self.leaf1, self.mid1),
            call(self.mid1, self.root),
            call(self.root, None),
        ])
