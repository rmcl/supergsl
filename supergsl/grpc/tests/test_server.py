from mock import patch, Mock
import unittest
from supergsl.grpc.server import SuperGSLCompilerService
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
    CreateCompilerSessionResult,
)


class SuperGSLCompilerServiceTestCase(unittest.TestCase):
    """Test case for `SuperGSLCompilerService` gRPC server."""

    def setUp(self):
        self.settings = {
            'YO': 'DUDE!'
        }
        self.service = SuperGSLCompilerService(self.settings)

    def test_CreateCompilerSession_creates_new_session(self):
        """Test that CreateCompilerSession correctly initializes a CompilerPipeline."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_mock:
            compiler_pipeline_mock.return_value = 'HELLO WORLD!'
            result = self.service.CreateCompilerSession(
                CreateCompilerSessionRequest(), Mock())

        compiler_pipeline_mock.assert_called_once_with(self.settings)

        session = self.service.get_compiler_session(result.session_identifier)
        self.assertEqual(session, 'HELLO WORLD!')
