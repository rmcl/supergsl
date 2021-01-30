from .lexer import Lexer
from .parser import ParserBuilder

from supergsl.core.symbol_table import SymbolTable
from supergsl.core.function import InvokeFunctionPass
from supergsl.core.plugin import PluginProvider
from supergsl.core.imports import ResolveImportsPass
from supergsl.core.parts.slice import ResolvePartSlicePass


class CompilerPipeline(object):

    def __init__(self):
        self._symbol_table = SymbolTable()
        self.plugins = PluginProvider(self._symbol_table)

    def get_backend_passes(self):
        return [
            ResolveImportsPass,
            ResolvePartSlicePass,
            InvokeFunctionPass,
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
            backend_pass_inst = backend_pass_class(self._symbol_table)

            print('performing pass... %s' % backend_pass_inst.get_pass_name())
            ast = backend_pass_inst.perform(ast)

        return ast

    def compile(self, source_code):
        ast = self.perform_frontend_compile(source_code)
        return self.perform_backend_compile(ast)

    def get_lexer(self):
        return Lexer().get_lexer()

    def get_parser(self):
        return ParserBuilder()
