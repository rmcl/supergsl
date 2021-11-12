"""Tests for SuperGSL's builtin gRPC server."""
from unittest import TestCase
from unittest.mock import patch, Mock
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.grpc.server import SuperGSLCompilerService
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
    CreateCompilerSessionResult,

    GetSequenceDetailRequest,
    GetSequenceDetailResult,

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

        self.symbol_table = self.fixtures.mk_symbol_table()
        self.symbol_table.insert('sequences', self.fixtures.sequence_store)

    def test_CreateCompilerSession_creates_new_session(self):
        """Test that CreateCompilerSession correctly initializes a CompilerPipeline."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_mock:
            compiler_pipeline_mock.return_value = 'HELLO WORLD!'
            result = self.service.CreateCompilerSession(
                CreateCompilerSessionRequest(), Mock())

        compiler_pipeline_mock.assert_called_once_with(self.settings)

        session = self.service.get_compiler_session(result.session_identifier)
        self.assertEqual(session, 'HELLO WORLD!')

    def test_ListSymbolTable_returns_symbols(self):
        """Verify that we get expected symbols from the symbol table."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_class_mock:
            compiler_pipeline = compiler_pipeline_class_mock.return_value
            compiler_pipeline.symbols = self.symbol_table

            identifier, compiler_mock_inst = self.service.create_compiler_session()

            result = self.service.ListSymbolTable(
                ListSymbolTableRequest(session_identifier=identifier), Mock())

        compiler_pipeline_class_mock.assert_called_once_with(self.settings)

        self.assertEqual(len(result.symbols), 3) # Two parts + the Sequence entry
        self.assertEqual(result.symbols['uHO'].type, 'Part')
        self.assertEqual(result.symbols['awesome.tHUG'].type, 'Part')

    def test_GetSequenceDetail_include_sequence(self):
        """GetSequenceDetail including the sequence of the entry."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_class_mock:
            compiler_pipeline = compiler_pipeline_class_mock.return_value
            compiler_pipeline.symbols = self.symbol_table

            session_identifier, compiler_mock_inst = self.service.create_compiler_session()

            sequence_entry = self.fixtures.mk_random_dna_sequence_entry(150)
            sequence_entry.roles = [self.fixtures.mk_sequence_role()]

            result = self.service.GetSequenceDetail(
                GetSequenceDetailRequest(
                    session_identifier=session_identifier,
                    sequence_identifier=str(sequence_entry.id),
                    include_sequence=True),
                Mock())

            self.assertEqual(result.sequence_entry.sequence, str(sequence_entry.sequence))

    def test_GetSequenceDetail(self):
        """Retrieve details about a sequence in the Sequence Store."""
        with patch('supergsl.grpc.server.CompilerPipeline') as compiler_pipeline_class_mock:
            compiler_pipeline = compiler_pipeline_class_mock.return_value
            compiler_pipeline.symbols = self.symbol_table

            session_identifier, compiler_mock_inst = self.service.create_compiler_session()

            sequence_entry = self.fixtures.mk_random_dna_sequence_entry(150)
            sequence_entry.roles = [self.fixtures.mk_sequence_role()]

            result = self.service.GetSequenceDetail(
                GetSequenceDetailRequest(
                    session_identifier=session_identifier,
                    sequence_identifier=str(sequence_entry.id)),
                Mock())

        self.assertEqual(result.sequence_entry.identifier, str(sequence_entry.id))
        self.assertEqual(result.sequence_entry.is_composite, False)
        self.assertEqual([
            (role.uri, role.name, role.description)
            for role in result.sequence_entry.roles
        ], [
            ('test/role/1', 'TestRole', 'test role')
        ])
