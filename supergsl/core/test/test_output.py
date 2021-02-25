"""Tests for the output core module."""
import unittest
from mock import Mock
from io import StringIO
from .fixtures import SuperGSLCoreFixtures
from supergsl.core.output import (
    ASTPrintOutputProvider,
    OutputPipeline,
    TestOutputProvider
)


class OutputTestCase(unittest.TestCase):
    """Test that the output providers."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.example_ast = self.fixtures.get_assembly_ast()

    def test_output_resolve_providers(self):
        """Test that the OutputProvider resolve providers provided in the settings."""
        pipeline = OutputPipeline({
            'output_providers': [
                'supergsl.core.output.TestOutputProvider'
            ]
        })

        self.assertEqual(
            pipeline.get_available_outputers(), {
                'test': TestOutputProvider
            })

    def test_run_providers(self):
        pipeline = OutputPipeline({
            'output_providers': [
                'supergsl.core.output.TestOutputProvider'
            ]
        })

        one_inst = Mock()
        two_inst = Mock()
        pipeline.available_outputers = {
            'one': Mock(return_value=one_inst),
            'two': Mock(return_value=two_inst)
        }

        args = Mock()
        args.output_format = ['one', 'two']
        pipeline.validate_args(args)

        ast_mock = Mock()
        pipeline.run(ast_mock, args)

        one_inst.perform.assert_called_once_with(ast_mock)
        two_inst.perform.assert_called_once_with(ast_mock)

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
