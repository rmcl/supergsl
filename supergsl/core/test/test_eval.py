"""Tests for the eval module."""
import unittest
from mock import Mock
from supergsl.core.eval import DependencyGraph

class DependencyGraphTestCase(unittest.TestCase):
    """Test the DependencyGraph class."""

    def setUp(self):
        pass

    def test_dependency_graph(self):
        """Build up a dependency graph and then incremently remove nodes."""
        graph = DependencyGraph()

        nodes = [
            Mock(name=str(idx))
            for idx in range(5)
        ]

        graph.add_edge(nodes[0], nodes[1])
        graph.add_edge(nodes[0], nodes[2])
        graph.add_edge(nodes[0], nodes[3])

        graph.add_edge(nodes[1], nodes[4])

        self.assertEqual(graph.get_available_nodes(), [
            nodes[2],
            nodes[3],
            nodes[4]
        ])

        graph.remove_node(nodes[4])

        self.assertEqual(graph.get_available_nodes(), [
            nodes[1],
            nodes[2],
            nodes[3]
        ])

        graph.remove_node(nodes[1])
        graph.remove_node(nodes[2])

        self.assertEqual(graph.get_available_nodes(), [
            nodes[3]
        ])

        graph.remove_node(nodes[3])

        self.assertEqual(graph.get_available_nodes(), [
            nodes[0]
        ])

        graph.remove_node(nodes[0])

        self.assertEqual(graph.get_available_nodes(), [])

    def test_dependency_graph_get_available_node(self):
        """Add several nodes and retrieve those that lack dependency edges."""
        graph = DependencyGraph()

        nodes = [
            Mock(name=idx)
            for idx in range(10)
        ]

        for idx in range(1, len(nodes)):
            node = nodes[idx-1]
            dep_node = nodes[idx]

            graph.add_edge(node, dep_node)

        last_node = nodes[-1]
        self.assertEqual(graph.get_available_nodes(), [last_node])

        graph.remove_node(last_node)

        self.assertEqual(graph.get_available_nodes(), [nodes[-2]])
