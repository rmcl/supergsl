from rply import LexerGenerator


class Lexer():
    def __init__(self):
        self.lexer = LexerGenerator()

    def _add_tokens(self):

        # Reserved Words
        self.lexer.add('FROM', r'from')
        self.lexer.add('IMPORT', r'import')

        # Other characters
        self.lexer.add('OPEN_PAREN', r'\(')
        self.lexer.add('CLOSE_PAREN', r'\)')
        self.lexer.add('OPEN_CURLY_BRACKET', r'\{')
        self.lexer.add('CLOSE_CURLY_BRACKET', r'\}')
        self.lexer.add('OPEN_BRACKET', r'\[')
        self.lexer.add('CLOSE_BRACKET', r'\]')
        self.lexer.add('COLON', r'\:')
        self.lexer.add('SEMICOLON', r'\;')
        self.lexer.add('COMMA', r'\,')
        self.lexer.add('PERIOD', r'\.')
        #self.lexer.add('HASH', r'\#')

        self.lexer.add('NUMBER', r'\d+')
        self.lexer.add('IDENTIFIER', r'\w[\w\d\_\-]*')

        # Ignore spaces
        self.lexer.ignore('\s+')

    def get_lexer(self):
        self._add_tokens()
        return self.lexer.build()
