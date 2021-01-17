from graphviz import Digraph

from supergsl.core.backend import BreadthFirstNodeFilteredPass


class ASTGraphPass(BreadthFirstNodeFilteredPass):
    name = 'ASTDottyGraph'
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
        self.ast_graph = Digraph(comment='AST #1')
        return ast

    def visit(self, node, parent_node):
        node_name = self.get_node_name(node)
        print('NODE', node, 'hi', node_name)
        self.ast_graph.node(node_name, node.get_node_label())

        parent_node_name = self.get_node_name(parent_node)
        self.ast_graph.edge(parent_node_name, node_name)


    def after_pass(self, ast):
        print(self.ast_graph.source)
        return ast


    '''
    def program_import_list(self, node):
        #imports = Graph(name='Program Imports', node_attr={'shape': 'box'})
        #c.edge('foo', 'bar')


    def program_handler(self, ast_node):
        dot_code = """
            digraph G {
               subgraph imports {
                   style=filled;
                   color=lightgrey;
                   node [style=filled,color=white];
        """



        dot_code += """
                   label = "Imports";
               }
            }
        """

        return dot_code


    def assembly_list_handler(self):
        return self.handle_node(ast_node.assemblies[0])



    def program_import_list_handler(self, ast_node):
        result = [
            'start [shape=Mdiamond];'
        ] + [
            'start -> %d;' % i
            for i, _ in enumerate(ast_node.program_imports)
        ]
        return '\n'.join(result)
    '''
