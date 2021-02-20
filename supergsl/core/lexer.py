from rply import LexerGenerator


class Lexer():
    def __init__(self):
        self.lexer = LexerGenerator()

    def _add_tokens(self):
        # Reserved Words
        self.lexer.add('FROM', r'from')
        self.lexer.add('IMPORT', r'import')
        self.lexer.add('AS', r'as')
        self.lexer.add('LET', r'let')

        # Other characters
        self.lexer.add('OPEN_PAREN', r'\(')
        self.lexer.add('CLOSE_PAREN', r'\)')
        self.lexer.add('OPEN_CURLY_BRACKET', r'\{')
        self.lexer.add('CLOSE_CURLY_BRACKET', r'\}')
        self.lexer.add('OPEN_BRACKET', r'\[')
        self.lexer.add('CLOSE_BRACKET', r'\]')
        self.lexer.add('FORWARD_SLASH', r'\/')
        #self.lexer.add('BACKWARD_SLASH', r'\\')
        self.lexer.add('COLON', r'\:')
        self.lexer.add('SEMICOLON', r'\;')
        self.lexer.add('COMMA', r'\,')
        self.lexer.add('PERIOD', r'\.')
        self.lexer.add('DOLLAR_SIGN', r'$')
        #self.lexer.add('HASH', r'\#')

        self.lexer.add('EQUAL', r'=')
        self.lexer.add('TILDE', r'~')
        self.lexer.add('EXCLAMATION', r'!')
        self.lexer.add('NUMBER', r'-?\d+')
        self.lexer.add('IDENTIFIER', r'\w[\w\d\_\-]*')

        # Ignore spaces
        self.lexer.ignore('\s+')

        # Comments - Ignore multiline comments in c syntax
        # Example: /* This is a comment! */
        self.lexer.ignore(r'/\*([\s\S]*?)\*/\s*')

        # Comments - Ignore remainder of line starting with "#".
        self.lexer.ignore(r'#.*\n')

    def get_lexer(self):
        self._add_tokens()
        return self.lexer.build()
