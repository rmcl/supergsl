from graphviz import Digraph, Graph

from supergsl.core.backend import BreadthFirstNodeFilteredPass


class PartSliceTreePass(BreadthFirstNodeFilteredPass):

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly,
        }

    def get_node_name(self, node):
        try:
            return self.node_names[node]
        except KeyError:
            n = self._node_names[self.cur_node_count]
            self.cur_node_count += 1
            self.node_names[node] = n
            return n

    def before_pass(self, ast):
        self.node_names = {}
        self._node_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.cur_node_count = 0
        self.graph = Digraph(comment='Parts')
        return ast

    def after_pass(self, ast):
        print(self.graph.source)
        return ast

    def visit_part(self, part_graph, part):
        part_node_name = self.get_node_name(part)
        part_graph.node(part_node_name, part.identifier)

        if part.parent_part:
            self.visit_part(part_graph, part.parent_part)
            parent_part_node_name = self.get_node_name(part.parent_part)
            part_graph.edge(part_node_name, parent_part_node_name)
        elif part.start.reference:
            reference = part.provider
            reference_node_name = self.get_node_name(reference)
            part_graph.node(reference_node_name, reference.name)
            part_graph.edge(part_node_name, reference_node_name)

    def visit_assembly(self, assembly_node):
        assembly_graph = Digraph(name=str(assembly_node))

        for part_node in assembly_node.parts:
            part_graph = Digraph(name=part_node.identifier)
            part = part_node.part

            self.visit_part(part_graph, part)

            assembly_graph.subgraph(part_graph)

        self.graph.subgraph(assembly_graph)
        return assembly_node

class ASTGraphPass(BreadthFirstNodeFilteredPass):
    name = 'ASTDotGraph'
    allow_modification = False

    def get_node_name(self, node):
        try:
            return self.node_names[node]
        except KeyError:
            n = self._node_names[self.cur_node_count]
            self.cur_node_count += 1
            self.node_names[node] = n
            return n

    def before_pass(self, ast):
        self.node_names = {}
        self._node_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.cur_node_count = 0
        self.ast_graph = Digraph(comment='AST #1', node_attr={'shape': 'box'})
        return ast

    def visit(self, node, parent_node):
        node_name = self.get_node_name(node)
        parent_node_name = self.get_node_name(parent_node)

        self.ast_graph.node(node_name, node.get_node_label())
        self.ast_graph.edge(parent_node_name, node_name)

    def after_pass(self, ast):
        print(self.ast_graph.source)
        return ast
