"""Unit tests for the Compiler Pipeline."""
from unittest import TestCase
from unittest.mock import Mock, patch, call
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

    def test_register_provider(self):
        """Validate that we can register a provider from a pipeline instance."""
        pipeline = CompilerPipeline({
            'plugins': []
        })

        provider_class = Mock()
        provider_inst = pipeline.register_provider(
            'this.path', provider_class, arg1='hello', arg2='hi')

        self.assertEqual(provider_class.return_value, provider_inst)
        provider_class.assert_called_once()
