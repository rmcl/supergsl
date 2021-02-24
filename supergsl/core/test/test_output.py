"""Tests for the output core module."""
import unittest
from io import StringIO
from .fixtures import SuperGSLCoreFixtures
from supergsl.core.output import (
    ASTPrintOutputProvider
)


class OutputTestCase(unittest.TestCase):
    """Test that the output providers."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.example_ast = self.fixtures.get_assembly_ast()

    def test_primer_output_provider(self):
        """Test that ASTPrintOutput prints the AST as expected."""
        output_provider = ASTPrintOutputProvider(None, allow_modification=False)

        output_provider.stream = StringIO()
        output_provider.perform(self.example_ast)
        result = output_provider.stream.getvalue()

        self.assertEqual(
            result.replace(' ','').replace('\n', ''),
            str(self.example_ast_1).replace(' ', '').replace('\n', ''))


    example_ast_1 = {
        'label': None,
        'node': 'Assembly',
        'parts': [{'identifier': 'uHO',
            'invert': False,
            'node': 'Part',
            'slice': None},
           {'identifier': 'pADH1',
            'invert': False,
            'node': 'Part',
            'slice': None},
           {'identifier': 'gERG10',
            'invert': False,
            'node': 'Part',
            'slice': None},
           {'identifier': 'tADH1',
            'invert': False,
            'node': 'Part',
            'slice': None},
           {'identifier': 'dHO',
            'invert': False,
            'node': 'Part',
            'slice': None}]
    }
