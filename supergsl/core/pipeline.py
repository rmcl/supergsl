"""Define the core pipeline for compiling and running SuperGSL."""
from typing import Dict, List, Union, cast
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.types.builtin import SuperGSLType
from supergsl.core.plugin import PluginProvider
from supergsl.core.backend import BackendPipelinePass
from supergsl.core.eval import EvaluatePass
from supergsl.utils import resolve_import

from .lexer import SuperGSLLexer
from .parser import SuperGSLParser


class CompilerPipeline:
    """Orchestrate the conversion of superGSL source code to compiled sequences."""

    def __init__(self, settings):
        self._global_symbol_table = SymbolTable('global', None)
        self._settings = settings
        self.plugins = PluginProvider(self._global_symbol_table, self._settings)

    def compile(self, source_code):
        """Run the compiler on the provided source code."""
        ast = self.perform_frontend_compile(source_code)
        return self.perform_backend_compile(ast)


    def import_symbols(
        self,
        module_path : str,
        import_identifier_list : Union[str, List[str]]
    ) -> Union[
        SuperGSLType,
        Dict[str, SuperGSLType]
    ]:
        """Import a symbol from a provider into the SuperGSL symbol table.

        This is equivalent to the SuperGSL statement:
        ```
        from <provider> import a, b, c
        ```

        where provider is the module_path and
        import_identifier_list is a list of strings of format ['a', 'b', 'c']
        """

        if isinstance(import_identifier_list, str):
            import_identifier_list = [import_identifier_list]

        imported_symbols = {}
        for import_identifier in import_identifier_list:
            new_symbols = resolve_import(
                self._global_symbol_table,
                module_path.split('.'),
                import_identifier,
                None)

            imported_symbols.update(new_symbols)

        return imported_symbols


    @property
    def symbols(self) -> SymbolTable:
        """Access the compiler's global symbol table."""
        return self._global_symbol_table


    def get_backend_passes(self) -> List[BackendPipelinePass]:
        """Return an ordered list of compiler backend passes to be executed."""
        return cast(List[BackendPipelinePass], [
            EvaluatePass
        ])


    def get_provider(self, module_path):
        """Return the provider instantiated at a give module path."""
        import_table = self._global_symbol_table.enter_nested_scope('imports')
        return import_table.lookup(module_path)


    def perform_frontend_compile(self, source_code):
        """Generate an IR from SuperGSL source code.

        The compiler front end analyzes the source code to build an internal representation (IR)
        of the program which in the case of superGSL is a Abstract Syntax Tree (AST).

        Input: (str) source code
        Output: `ast.Program`
        """
        lexer = self.get_lexer()
        tokens = lexer.lex(source_code)

        parser = self.get_parser()
        return parser.parse(tokens)


    def perform_backend_compile(self, ast):
        """Execute a series of backend compiler passes over the given AST."""
        pass_classes = self.get_backend_passes()

        for backend_pass_class in pass_classes:
            backend_pass_inst = backend_pass_class(self._global_symbol_table)

            ast = backend_pass_inst.perform(ast)

        return ast


    def get_lexer(self):
        """Retrieve a reference to the lexer."""
        return SuperGSLLexer()


    def get_parser(self):
        """Retrieve a reference to the parser."""
        return SuperGSLParser.create_supergsl_parser()
