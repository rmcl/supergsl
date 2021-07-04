"""Unit tests for the Compiler Pipeline."""
from unittest import TestCase
from unittest.mock import patch, call
from supergsl.core.pipeline import CompilerPipeline


class CompilerPipelineTestCase(TestCase):
    """Test the compiler pipeline."""
    maxDiff = None

    @patch('supergsl.core.pipeline.resolve_import')
    def test_import_symbols(self, resolve_import_mock):
        """The compiler pipeline can import symbols."""
        resolve_import_mock.return_value = {
            'HELLO': 'world'
        }

        pipeline = CompilerPipeline({
            'plugins': {}
        })

        results = pipeline.import_symbols('test.test', [
            'a', 'b', 'c'
        ])

        self.assertEqual(results, {'HELLO': 'world'})
        resolve_import_mock.assert_has_calls([
            call(pipeline.symbols, ['test', 'test'], 'a', None),
            call(pipeline.symbols, ['test', 'test'], 'b', None),
            call(pipeline.symbols, ['test', 'test'], 'c', None)
        ])
