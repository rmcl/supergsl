"""Unit tests for the SuperGSL parser."""
import unittest
import mock
from rply import Token
from supergsl.core.parser import ParserBuilder
from supergsl.core.ast import Program
from supergsl.core.exception import ParsingError


class ParserTestCase(unittest.TestCase):
    """Test that the parser correctly tokens into a valid AST."""
    maxDiff = None

    def setUp(self):
        self.parser = ParserBuilder()

    def test_build_ast_import(self):
        """Test building an AST from the parsed tokens of "from S288C import ADHA, ERG10, HO"."""

        tokens = iter((
            # Import
            Token('FROM', 'from'),
            Token('IDENTIFIER', 'S288C'),
            Token('IMPORT', 'import'),
            Token('IDENTIFIER', 'ADHA'),
            Token('AS', 'as'),
            Token('IDENTIFIER', 'BLOOP'),
            Token('COMMA', ','),
            Token('IDENTIFIER', 'ERG10'),
            Token('COMMA', ','),
            Token('IDENTIFIER', 'HO'),

            # Single Assembly
            Token('IDENTIFIER', 'uHO'),
        ))
        ast = self.parser.parse(tokens)

        self.assertEquals(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'items': [{
                    'label': None,
                    'node': 'Assembly',
                    'parts': [{
                        'node': 'Part',
                        'identifier': 'uHO',
                        'invert': False,
                        'slice': None
                    }]
                }],
                'node': 'DefinitionList'
            },
            'imports': [{
                'imports': [
                    {'alias': 'BLOOP', 'identifier': 'ADHA'},
                    {'alias': None, 'identifier': 'ERG10'},
                    {'alias': None, 'identifier': 'HO'}
                ],
                'module': ['S288C'],
                'node': 'Import'
            }],
            'node': 'Program'
        })

    def test_build_ast_assembly(self):
        tokens = iter((
            Token('IDENTIFIER', 'uHO'),
            Token('SEMICOLON', ';'),
            Token('IDENTIFIER', 'pADH1')
        ))
        ast = self.parser.parse(tokens)

        self.assertEqual(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'node': 'DefinitionList',
                'items': [{
                    'label': None,
                    'node': 'Assembly',
                    'parts': [{
                        'node': 'Part',
                        'identifier': 'uHO',
                        'invert': False,
                        'slice': None
                    }, {
                        'node': 'Part',
                        'identifier': 'pADH1',
                        'invert': False,
                        'slice': None
                    }],
                }]
            },
            'imports': [],
            'node': 'Program'
        })

    def test_build_ast_function_call(self):
        tokens = iter((
            Token('IDENTIFIER', 'functest'),
            Token('OPEN_CURLY_BRACKET', '{'),
            Token('IDENTIFIER', 'uHO'),
            Token('CLOSE_CURLY_BRACKET', '}')
        ))
        ast = self.parser.parse(tokens)

        self.assertEqual(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'node': 'DefinitionList',
                'items': [{
                    'identifier': 'functest',
                    'node': 'FunctionInvocation',
                    'children': {
                        'node': 'DefinitionList',
                        'items': [{
                            'label': None,
                            'node': 'Assembly',
                            'parts': [{
                                'node': 'Part',
                                'identifier': 'uHO',
                                'invert': False,
                                'slice': None
                            }]
                        }]
                    },
                    'params': None,
                    'label': None
                }]
            },
            'imports': [],
            'node': 'Program'
        })

    def test_build_ast_function_call_with_empty_params(self):
        """Parse a AST containing a function with a nested body into an AST."""
        tokens = iter((
            Token('IDENTIFIER', 'functest'),
            Token('OPEN_CURLY_BRACKET', '{'),
            Token('IDENTIFIER', 'uHO'),
            Token('CLOSE_CURLY_BRACKET', '}')
        ))
        ast = self.parser.parse(tokens)

        self.assertEqual(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'node': 'DefinitionList',
                'items': [{
                    'identifier': 'functest',
                    'node': 'FunctionInvocation',
                    'children': {
                        'node': 'DefinitionList',
                        'items': [{
                            'label': None,
                            'node': 'Assembly',
                            'parts': [{
                                'node': 'Part',
                                'identifier': 'uHO',
                                'invert': False,
                                'slice': None
                            }]
                        }]
                    },
                    'params': None,
                    'label': None
                }]
            },
            'imports': [],
            'node': 'Program'
        })

    def test_build_ast_function_call_with_params(self):
        """Parse a function call into a AST."""
        tokens = iter((
            Token('IDENTIFIER', 'functest'),
            Token('OPEN_PAREN', '('),
            Token('IDENTIFIER', 'CHEESE'),
            Token('CLOSE_PAREN', ')'),
        ))
        ast = self.parser.parse(tokens)

        self.assertEqual(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'node': 'DefinitionList',
                'items': [{
                    'identifier': 'functest',
                    'node': 'FunctionInvocation',
                    'children': None,
                    'params': [
                        'CHEESE'
                    ],
                    'label': None
                }]
            },
            'imports': [],
            'node': 'Program'
        })

    def test_build_ast_sequence_constant(self):
        """Parse to an AST for a constant nucleotide sequence."""
        tokens = iter((
            Token('FORWARD_SLASH', '/'),
            Token('IDENTIFIER', 'ATGG'),
            Token('FORWARD_SLASH', '/'),

        ))
        ast = self.parser.parse(tokens)

        self.assertEqual(type(ast), Program)
        self.assertEqual(ast.eval(), {
            'definitions': {
                'node': 'DefinitionList',
                'items': [{
                    'label': None,
                    'node': 'Assembly',
                    'parts': [{
                        'node': 'SequenceConstant',
                        'sequence': 'ATGG',
                        'type': 'DNA'
                    }]
                }],
            },
            'imports': [],
            'node': 'Program'
        })
