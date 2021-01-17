from rply.errors import LexingError

from .lexer import Lexer
from .parser import ParserBuilder

from supergsl.core.ast import SymbolRepository
from supergsl.core.parts import PartSymbolTable
from supergsl.core.function import FunctionSymbolTable, InvokeFunctionPass
from supergsl.core.plugin import PluginProvider
from supergsl.core.imports import ResolveImportsPass
from supergsl.core.parts.slice import ResolvePartSlicePass

from supergsl.plugins.graphviz import ASTGraphPass

class CompilerPipeline(object):

    def __init__(self):
        self.initialize_symbol_repository()

        self.plugins = PluginProvider(self._symbol_registry)

    def initialize_symbol_repository(self):
        self._symbol_registry = SymbolRepository()
        self._symbol_registry.register('parts', PartSymbolTable())
        self._symbol_registry.register('functions', FunctionSymbolTable())

    def get_backend_passes(self):
        return [
            ResolveImportsPass,
            ResolvePartSlicePass,
            ASTGraphPass,
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

        try:
            return parser.parse(tokens)
        except LexingError as error:

            ## Todo(rmcl): BIG TODO HERE! How do we want to display lex errors?
            error_pos = error.getsourcepos().idx
            context_start = min(error_pos-25, 0)
            context_end = max(error_pos+50, len(source_code))
            code_context = '%s**%s**%s' % (
                source_code[context_start:error_pos],
                source_code[error_pos],
                source_code[error_pos+1:context_end]
            )


            print('Syntax error: %s' % code_context)

    def perform_backend_compile(self, ast):
        pass_classes = self.get_backend_passes()

        print('BACKEND!!!')

        for backend_pass_class in pass_classes:
            backend_pass_inst = backend_pass_class(self._symbol_registry)

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
