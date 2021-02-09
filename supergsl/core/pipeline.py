from typing import List, cast
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.function import InvokeFunctionPass
from supergsl.core.plugin import PluginProvider
from supergsl.core.imports import ResolveImportsPass
from supergsl.core.backend import BackendPipelinePass
from supergsl.core.parts.slice import ResolvePartSlicePass

from .lexer import Lexer
from .parser import ParserBuilder


class CompilerPipeline(object):
    """Orchestrate the conversion of superGSL source code to compiled sequences."""

    def __init__(self, settings):
        self._symbol_table = SymbolTable()
        self._settings = settings
        self.plugins = PluginProvider(self._symbol_table, self._settings)

    def get_backend_passes(self) -> List[BackendPipelinePass]:
        """Return an ordered list of compiler backend passes to be executed."""
        return cast(List[BackendPipelinePass], [
            ResolveImportsPass,
            ResolvePartSlicePass,
            InvokeFunctionPass,
        ])

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
        """Execute a series of backend compiler passes over the given AST."""
        pass_classes = self.get_backend_passes()

        for backend_pass_class in pass_classes:
            backend_pass_inst = backend_pass_class(self._symbol_table)

            print('performing pass... %s' % backend_pass_inst.get_pass_name())
            ast = backend_pass_inst.perform(ast)

        return ast

    def compile(self, source_code):
        """Run the compiler on the provided source code."""
        ast = self.perform_frontend_compile(source_code)
        return self.perform_backend_compile(ast)

    def get_lexer(self):
        """Retrieve a reference to the lexer."""
        return Lexer().get_lexer()

    def get_parser(self):
        """Retrieve a reference to the parser."""
        return ParserBuilder()
