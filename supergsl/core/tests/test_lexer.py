import unittest
from rply.errors import LexingError
from supergsl.core.lexer import SuperGSLLexer


class LexerTestCase(unittest.TestCase):
    """Test that the lexer correctly translates source code into tokens."""
    maxDiff = None

    def setUp(self):
        self.lexer = SuperGSLLexer()

    def test_single_quote_multiline_string_constant_not_allowed(self):
        """Newline characters are not allowed in single quoted constants."""
        input_string = """
            'HELLO
            WORLD'
        """
        with self.assertRaises(LexingError):
            list(self.lexer.lex(input_string))

    def test_string_constant_does_not_match_bad_constants(self):
        """"""
        bad_examples = [
            "\'\'\'",
        ]

        for example_str in bad_examples:
            with self.assertRaises(LexingError):
                list(self.lexer.lex(example_str))


    def test_scan_some_examples(self):
        """Test that lexer output matches expectation for a number of examples."""
        for idx, example in enumerate(examples):
            inp = example[0]
            expected_tokens = list(example[1:])

            try:
                tokens = [
                    (
                        token.gettokentype(),
                        token.value
                    )
                    for token in self.lexer.lex(inp)
                ]
            except LexingError as error:
                print(inp)
                raise error


            self.assertEqual(
                tokens,
                expected_tokens,
                'Unexpected tokens for example %d' % idx)

examples = [
    (
        "uHO ; pADH1",
        ('IDENTIFIER', 'uHO'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'pADH1')

    ), (
        "uHO ; pADH1 ; gERG10[1:728] ; dHO",
        ('IDENTIFIER', 'uHO'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'pADH1'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'gERG10'),
        ('OPEN_BRACKET', '['),
        ('NUMBER', '1'),
        ('COLON', ':'),
        ('NUMBER', '728'),
        ('CLOSE_BRACKET', ']'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'dHO')
    ), (
        """
        from S288C import ADHA as ADHA_is_great, ERG10, HO
        uHO ; pADH1 ; gERG10[1:728] ; dHO
        """,
        ('FROM', 'from'),
        ('IDENTIFIER', 'S288C'),
        ('IMPORT', 'import'),
        ('IDENTIFIER', 'ADHA'),
        ('AS', 'as '),
        ('IDENTIFIER', 'ADHA_is_great'),
        ('COMMA', ','),
        ('IDENTIFIER', 'ERG10'),
        ('COMMA', ','),
        ('IDENTIFIER', 'HO'),

        ('IDENTIFIER', 'uHO'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'pADH1'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'gERG10'),
        ('OPEN_BRACKET', '['),
        ('NUMBER', '1'),
        ('COLON', ':'),
        ('NUMBER', '728'),
        ('CLOSE_BRACKET', ']'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'dHO')
    ), (
        "/* HELLO THIS IS A \n\n COMMENT */ uHO ; pADH1",
        ('IDENTIFIER', 'uHO'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'pADH1')
    ), (
        "/* COMMENT */ uHO ; /* COMMENT \n\n\n COMMENT */ pADH1 /* COMMENT */",
        ('IDENTIFIER', 'uHO'),
        ('SEMICOLON', ';'),
        ('IDENTIFIER', 'pADH1')
    ), (
        "/$MAAADR*/",
        ('FORWARD_SLASH', '/'),
        ('AMINO_ACID_SEQUENCE', '$MAAADR*'),
        ('FORWARD_SLASH', '/'),
    ), (
        "'HELLO WORLD'",
        ('STRING_CONSTANT', '\'HELLO WORLD\'')
    ), (
        "'HELLO' 'HELO'",
        ('STRING_CONSTANT', '\'HELLO\''),
        ('STRING_CONSTANT', '\'HELO\'')
    )
]
