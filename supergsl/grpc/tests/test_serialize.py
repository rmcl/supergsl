from unittest import TestCase
from unittest.mock import patch, Mock
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
)
from supergsl.grpc.serialize import role_dto


class SuperGSLCompilerServiceTestCase(TestCase):
    """Test case for `SuperGSLCompilerService` gRPC server."""

    def setUp(self):
        self.settings = {
            'YO': 'DUDE!'
        }
        self.fixtures = SuperGSLCoreFixtures()

        self.symbol_table = self.fixtures.mk_symbol_table()
        self.symbol_table.insert('sequences', self.fixtures.sequence_store)

    def test_role_dto(self):
        """Test the role dto function."""
        role = self.fixtures.mk_sequence_role()
        role_grpc_obj = role_dto(role)

        self.assertEqual(role_grpc_obj.uri, role.uri)
        self.assertEqual(role_grpc_obj.name, role.name)
        self.assertEqual(role_grpc_obj.description, role.description)
