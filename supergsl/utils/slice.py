from supergsl.core.lexer import SuperGSLLexer
from supergsl.core.parser import SuperGSLParser
from supergsl.core.eval import EvaluatePass


def parse_slice_str(slice_source_code : str) -> 'Slice':
    lexer = SuperGSLLexer()
    tokens = lexer.lex(slice_source_code)

    parser = SuperGSLParser.create_slice_parser()
    ast = parser.parse(tokens)

    eval_pass = EvaluatePass(None)
    return eval_pass.visit(ast)
