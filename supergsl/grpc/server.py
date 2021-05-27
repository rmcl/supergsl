import uuid
from typing import Optional, Dict
from concurrent import futures
from grpc import server as grpc_server

from supergsl.grpc.stubs.sgsl_pb2_grpc import (
    SuperGSLCompilerServicer,
    add_SuperGSLCompilerServicer_to_server
)
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
    CreateCompilerSessionResult,
    ListSymbolTableRequest,
    ListSymbolTableResult,
    CompileRequest,
    CompileResult,
    Symbol
)

from supergsl.core.pipeline import CompilerPipeline


class SuperGSLCompilerService(SuperGSLCompilerServicer):
    """Implement the SuperGSL Compiler gRPC Service."""

    def __init__(self, compiler_settings : dict):
        self.compiler_settings = compiler_settings
        self.server : Optional[grpc_server] = None

        self.compiler_sessions : Dict[str, CompilerPipeline] = {}

    def CreateCompilerSession(
        self,
        request : CreateCompilerSessionRequest,
        context
    ) -> CreateCompilerSessionResult:
        """Missing associated documentation comment in .proto file."""

        identifier = str(uuid.uuid4())
        self.compiler_sessions[identifier] = CompilerPipeline(self.compiler_settings)
        return CreateCompilerSessionResult(session_identifier=identifier)

    def ListSymbolTable(
        self, request : ListSymbolTableRequest,
        context
    ) -> ListSymbolTableResult:
        """List the symbols from the given compiler session's symbol table."""
        print(type(request))
        print(type(context))
        compiler_pipeline = self.get_compiler_session(request.session_identifier)

        symbol_table = compiler_pipeline.get_symbol_table()
        symbol_map = {}
        for key, value in symbol_table._symbols.items():
            symbol = Symbol(type=type(value).__name__)
            symbol_map[key] = symbol

        return ListSymbolTableResult(symbols=symbol_map)

    def Compile(
        self,
        request : CompileRequest,
        context
    ) -> CompileResult:
        """Compile a string of source code in the specified compiler session."""
        compiler_pipeline = self.get_compiler_session(request.session_identifier)
        compiler_pipeline.compile(request.source_code)
        return CompileResult(
            success=True,
            error='test'
        )

    def get_compiler_session(self, session_identifier : str) -> CompilerPipeline:
        """Return a CompilerPipeline corresponding to the provided session identifier."""
        return self.compiler_sessions[session_identifier]

    def start_listening(self):
        """Start the gRPC server in a new thread."""
        self.server = grpc_server(futures.ThreadPoolExecutor(max_workers=1))
        add_SuperGSLCompilerServicer_to_server(self, self.server)
        self.server.add_insecure_port('[::]:50051')
        self.server.start()
        self.server.wait_for_termination()
