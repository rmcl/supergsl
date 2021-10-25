import uuid
from typing import Optional, Dict, Tuple, Any
from concurrent import futures
from grpc import server as grpc_server
from google.protobuf.struct_pb2 import Struct

from supergsl.grpc.stubs.sgsl_pb2_grpc import (
    SuperGSLCompilerServicer,
    add_SuperGSLCompilerServicer_to_server
)
from supergsl.grpc.stubs.sgsl_pb2 import (
    CreateCompilerSessionRequest,
    CreateCompilerSessionResult,
    ListSequenceStoreRequest,
    ListSequenceStoreResult,
    ListSymbolTableRequest,
    ListSymbolTableResult,
    SequenceEntry,
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

        identifier, _ = self.create_compiler_session()
        return CreateCompilerSessionResult(session_identifier=identifier)

    def serialize_symbol_table(self, symbol_table, prefix=None) -> Dict[str, Any]:
        """Recursively convert symbol table items to "dot" notation."""
        result = {}
        for key, value in symbol_table:
            if prefix:
                full_key_path = '%s.%s' % (prefix, key)
            else:
                full_key_path = key

            result[full_key_path] = value

        for nested_scope_key, nested_scope_table in symbol_table.nested_scopes():
            if prefix:
                new_prefix = '%s.%s' % (prefix, nested_scope_key)
            else:
                new_prefix = nested_scope_key

            result.update(self.serialize_symbol_table(nested_scope_table, new_prefix))

        return result


    def ListSequenceStore(
        self,
        request : ListSequenceStoreRequest,
        context
    ):
        """Return a map of Sequence Entries in the compiler Sequence Store."""
        serialized_entries : Dict[str, SequenceEntry] = {}

        compiler_pipeline = self.get_compiler_session(request.session_identifier)
        sequence_store = compiler_pipeline.symbols['sequences']

        for identifier, sequence_entry in sequence_store.items():
            serialized_entry = SequenceEntry(
                identifier=str(identifier),
                is_composite=sequence_entry.is_composite,
                sequence=''
            )
            serialized_entries[str(identifier)] = serialized_entry

        return ListSequenceStoreResult(sequence_entries=serialized_entries)

    def GetSequenceDetail(self, request, context):
        """Missing associated documentation comment in .proto file."""
        compiler_pipeline = self.get_compiler_session(request.session_identifier)
        sequence_store = compiler_pipeline.symbols['sequences']
        sequence_entry = sequence_store.lookup(request.sequence_identifier)

        sequence = ''
        if request.include_sequence:
            sequence = sequence_entry.sequence

        entry = SequenceEntry(
            identifier=sequence_entry.id,
            is_composite=sequence_entry.is_composite,
            sequence=sequence
        )

        return GetSequenceDetailResult(sequence_entry=entry)

    def ListSymbolTable(
        self,
        request : ListSymbolTableRequest,
        context
    ) -> ListSymbolTableResult:
        """List the symbols from the given compiler session's symbol table."""
        compiler_pipeline = self.get_compiler_session(request.session_identifier)

        symbol_table = compiler_pipeline.symbols
        symbol_map : Dict[str, Symbol] = {}
        for key, value in self.serialize_symbol_table(symbol_table).items():
            details = Struct()
            details.update(value.serialize())

            symbol_map[key] = Symbol(
                type=type(value).__name__,
                details=details)

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

    def create_compiler_session(self) -> Tuple[str, CompilerPipeline]:
        """Create a Compiler session and store it."""
        identifier = str(uuid.uuid4())
        self.compiler_sessions[identifier] = CompilerPipeline(self.compiler_settings)
        return identifier, self.compiler_sessions[identifier]

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
