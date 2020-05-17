from supergsl.core.backend import BackendPipelinePass


class DottyASTGenerator(BackendPipelinePass):
    name = 'Dotty Generator Pass'

    def __init__(self):
        self.NODE_HANDLERS = {
            'Program': self.program_handler,
            'ProgramImportList': self.program_import_list_handler,

            'AssemblyList': self.assembly_list_handler,
        }

    def perform(self, ast):
        output = self.handle_node(ast)
        print(output)

    def handle_node(self, ast_node):
        node_type = type(ast_node).__name__
        handler = self.NODE_HANDLERS[node_type]

        return handler(ast_node)

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
