"""Define util functions for using supergsl lang slice strings."""
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.types.position import Slice

from supergsl.lang.lexer import SuperGSLLexer
from supergsl.lang.parser import SuperGSLParser
from supergsl.lang.eval import EvaluatePass



def parse_slice_str(slice_source_code : str) -> Slice:
    """Parse a part slice coordinate string into a `supergsl.core.types.slice.Slice`."""
    lexer = SuperGSLLexer()
    tokens = lexer.lex(slice_source_code)

    parser = SuperGSLParser.create_slice_parser()
    ast = parser.parse(tokens)

    eval_pass = EvaluatePass(SymbolTable('slice_eval', None))
    return eval_pass.visit(ast, SymbolTable('temp', None))
