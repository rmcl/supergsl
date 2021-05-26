import uuid
from typing import Optional
import grpc
from grpc import Server
from concurrent import futures
import supergsl.grpc.stubs.sgsl_pb2_grpc as pb2_grpc
import supergsl.grpc.stubs.sgsl_pb2 as pb2

from supergsl.core.pipeline import CompilerPipeline


class SuperGSLCompilerService(pb2_grpc.SuperGSLCompilerServicer):

    def __init__(self, compiler_settings):
        self.compiler_settings = compiler_settings
        self.server : Optional[Server] = None

        self.compiler_sessions = {}

    def CreateCompilerSession(self, request, context):
        """Missing associated documentation comment in .proto file."""

        identifier = str(uuid.uuid4())
        self.compiler_sessions[identifier] = CompilerPipeline(self.compiler_settings)
        return pb2.CreateCompilerSessionResult(session_identifier=identifier)

    def ListSymbolTable(self, request, context):
        compiler_pipeline = self.get_compiler_session(request.session_identifier)

        symbol_table = compiler_pipeline.get_symbol_table()
        symbol_map = {}
        for key, value in symbol_table._symbols.items():
            symbol = pb2.Symbol(type=type(value).__name__)
            symbol_map[key] = symbol

        return pb2.ListSymbolTableResult(symbols=symbol_map)

    def Compile(self, request, context):
        compiler_pipeline = self.get_compiler_session(request.session_identifier)
        compiler_pipeline.compile(request.source_code)
        return pb2.CompileResult(
            success=True,
            error='test'
        )

    def get_compiler_session(self, session_identifier):
        return self.compiler_sessions[session_identifier]

    def start_listening(self):
        """Start the gRPC server in a new thread."""
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=1))
        pb2_grpc.add_SuperGSLCompilerServicer_to_server(self, self.server)
        self.server.add_insecure_port('[::]:50051')
        self.server.start()
        self.server.wait_for_termination()
