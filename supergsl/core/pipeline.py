from .lexer import Lexer
from .parser import ParserBuilder

class CompilerPipeline(object):

    def get_backend_passes(self):
        return [
            # DottyASTGenerator
        ]

    def perform_frontend_compile(self, source_code):
        """Generate an IR from SuperGSL source code.

        The compiler front end analyzes the source code to build an internal representation (IR)
        of the program which in the case of superGSL is a Abstract Syntax Tree (AST).

        Input: (str) source code
        Output: `ast.Program`
        """
        tokens = self.get_lexer().lex(source_code)

        parser = self.get_parser()
        return parser.parse(tokens)

    def perform_backend_compile(self, ast):
        pass_classes = self.get_backend_passes()

        print('BACKEND!!!')

        for backend_pass_class in pass_classes:
            backend_pass_inst = backend_pass_class()

            print('performing pass... %s' % backend_pass_inst.name)
            ast = backend_pass_inst.perform(ast)

        return ast

    def compile(self, source_code):
        ast = self.perform_frontend_compile(source_code)

        return self.perform_backend_compile(ast).eval()

    def get_lexer(self):
        return Lexer().get_lexer()

    def get_parser(self):
        parser_generator = ParserBuilder()
        parser_generator.parse()
        return parser_generator.get_parser()
