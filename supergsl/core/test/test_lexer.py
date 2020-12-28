import unittest
from supergsl.core.lexer import Lexer


class LexerTestCase(unittest.TestCase):
    """Test that the lexer correctly translates source code into tokens."""
    maxDiff = None

    def setUp(self):
        self.lexer = Lexer().get_lexer()

    def test_scan_some_examples(self):
        for idx, example in enumerate(examples):
            inp = example[0]
            expected_tokens = list(example[1:])

        tokens = [
            (
                token.gettokentype(),
                token.value
            )
            for token in self.lexer.lex(inp)
        ]

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
        from S288C import ADHA, ERG10, HO
        uHO ; pADH1 ; gERG10[1:728] ; dHO
        """,
        ('FROM', 'from'),
        ('IDENTIFIER', 'S288C'),
        ('IMPORT', 'import'),
        ('IDENTIFIER', 'ADHA'),
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
    )
]
