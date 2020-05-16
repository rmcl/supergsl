from .lexer import Lexer
from .parser import ParserBuilder


class CompilerPipeline(object):

    def get_backend_passes(self):
        return [

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
        return parser.parse(tokens).eval()

    def perform_backend_compile(self, ast):
        passes = self.get_backend_passes()

        print('BACKEND!!!')

        return ast

    def compile(self, source_code):
        ast = self.perform_frontend_compile(source_code)

        return self.perform_backend_compile(ast)

    def get_lexer(self):
        return Lexer().get_lexer()

    def get_parser(self):
        parser_generator = ParserBuilder()
        parser_generator.parse()
        return parser_generator.get_parser()
