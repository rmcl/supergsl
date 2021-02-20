"""Tests for the output core module."""
import unittest
from csv import DictReader
from mock import Mock
from io import StringIO
from supergsl.core.output import (
    ASTPrintOutputProvider,
    PrimerOutputProvider
)
from .fixtures import SuperGSLCoreFixtures


class OutputTestCase(unittest.TestCase):
    """Test that the output providers."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.example_ast = self.fixtures.get_assembly_ast()

    def test_primer_output_provider(self):
        """Test that primer output creates CSV with primers."""
        output_provider = PrimerOutputProvider(None, allow_modification=False)

        class NoCloseStringIO(StringIO):
            """Subclass StringIO so we can inspect buffer before closing."""
            def close(self):
                pass

            def really_close(self):
                super().close()

        primer_file_fp = NoCloseStringIO()
        output_provider._open_primer_file = Mock(return_value=primer_file_fp)

        output_provider.perform(self.example_ast)

        csv_output_str = primer_file_fp.getvalue()
        primer_file_fp.close()

        reader = DictReader(StringIO(csv_output_str))
        for result in reader:
            primer_record = output_provider.primers[result['Part Identifier']]
            self.assertEqual(result['Forward Primer'], primer_record['Forward Primer'])
            self.assertEqual(result['Reverse Primer'], primer_record['Reverse Primer'])


    def test_ast_print_output_provider(self):
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
