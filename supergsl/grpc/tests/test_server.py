from unittest.mock import patch, Mock
from unittest import TestCase, skip

from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.grpc.server import SuperGSLCompilerService
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
    CreateCompilerSessionResult,

    ListSymbolTableRequest,
    ListSymbolTableResult
)


class SuperGSLCompilerServiceTestCase(TestCase):
    """Test case for `SuperGSLCompilerService` gRPC server."""

    def setUp(self):
        self.settings = {
            'YO': 'DUDE!'
        }
        self.service = SuperGSLCompilerService(self.settings)
        self.fixtures = SuperGSLCoreFixtures()

    def test_CreateCompilerSession_creates_new_session(self):
        """Test that CreateCompilerSession correctly initializes a CompilerPipeline."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_mock:
            compiler_pipeline_mock.return_value = 'HELLO WORLD!'
            result = self.service.CreateCompilerSession(
                CreateCompilerSessionRequest(), Mock())

        compiler_pipeline_mock.assert_called_once_with(self.settings)

        session = self.service.get_compiler_session(result.session_identifier)
        self.assertEqual(session, 'HELLO WORLD!')

    @skip('further debugging required.')
    def test_ListSymbolTable_returns_symbols(self):
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_class_mock:
            compiler_pipeline = compiler_pipeline_class_mock.return_value
            symbol_table = self.fixtures.mk_symbol_table()
            compiler_pipeline.symbols = symbol_table

            identifier, compiler_mock_inst = self.service.create_compiler_session()

            result = self.service.ListSymbolTable(
                ListSymbolTableRequest(session_identifier=identifier), Mock())

        compiler_pipeline_class_mock.assert_called_once_with(self.settings)

        self.assertEqual(len(result.symbols), 2)
        self.assertEqual(result.symbols['uHO'].type, 'Part')
        self.assertEqual(result.symbols['awesome.tHUG'].type, 'Part')
