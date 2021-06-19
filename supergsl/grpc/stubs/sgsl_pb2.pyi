"""
@generated by mypy-protobuf.  Do not edit manually!
isort:skip_file
"""
import abc
import builtins
import concurrent.futures
import google.protobuf.descriptor
import google.protobuf.internal.containers
import google.protobuf.message
import google.protobuf.service
import typing
import typing_extensions

DESCRIPTOR: google.protobuf.descriptor.FileDescriptor = ...

class CreateCompilerSessionRequest(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...

    def __init__(self,
        ) -> None: ...
global___CreateCompilerSessionRequest = CreateCompilerSessionRequest

class CreateCompilerSessionResult(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    SESSION_IDENTIFIER_FIELD_NUMBER: builtins.int
    session_identifier: typing.Text = ...

    def __init__(self,
        *,
        session_identifier : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"session_identifier",b"session_identifier"]) -> None: ...
global___CreateCompilerSessionResult = CreateCompilerSessionResult

class DestroyCompilerSessionRequest(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...

    def __init__(self,
        ) -> None: ...
global___DestroyCompilerSessionRequest = DestroyCompilerSessionRequest

class DestroyCompilerSessionResult(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    SUCCESS_FIELD_NUMBER: builtins.int
    ERROR_FIELD_NUMBER: builtins.int
    success: builtins.bool = ...
    error: typing.Text = ...

    def __init__(self,
        *,
        success : builtins.bool = ...,
        error : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"error",b"error",u"success",b"success"]) -> None: ...
global___DestroyCompilerSessionResult = DestroyCompilerSessionResult

class CompileRequest(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    SESSION_IDENTIFIER_FIELD_NUMBER: builtins.int
    SOURCE_CODE_FIELD_NUMBER: builtins.int
    session_identifier: typing.Text = ...
    source_code: typing.Text = ...

    def __init__(self,
        *,
        session_identifier : typing.Text = ...,
        source_code : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"session_identifier",b"session_identifier",u"source_code",b"source_code"]) -> None: ...
global___CompileRequest = CompileRequest

class CompileResult(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    SUCCESS_FIELD_NUMBER: builtins.int
    ERROR_FIELD_NUMBER: builtins.int
    success: builtins.bool = ...
    error: typing.Text = ...

    def __init__(self,
        *,
        success : builtins.bool = ...,
        error : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"error",b"error",u"success",b"success"]) -> None: ...
global___CompileResult = CompileResult

class ListSymbolTableRequest(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    SESSION_IDENTIFIER_FIELD_NUMBER: builtins.int
    session_identifier: typing.Text = ...

    def __init__(self,
        *,
        session_identifier : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"session_identifier",b"session_identifier"]) -> None: ...
global___ListSymbolTableRequest = ListSymbolTableRequest

class ListSymbolTableResult(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    class SymbolsEntry(google.protobuf.message.Message):
        DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
        KEY_FIELD_NUMBER: builtins.int
        VALUE_FIELD_NUMBER: builtins.int
        key: typing.Text = ...

        @property
        def value(self) -> global___Symbol: ...

        def __init__(self,
            *,
            key : typing.Text = ...,
            value : typing.Optional[global___Symbol] = ...,
            ) -> None: ...
        def HasField(self, field_name: typing_extensions.Literal[u"value",b"value"]) -> builtins.bool: ...
        def ClearField(self, field_name: typing_extensions.Literal[u"key",b"key",u"value",b"value"]) -> None: ...

    SYMBOLS_FIELD_NUMBER: builtins.int

    @property
    def symbols(self) -> google.protobuf.internal.containers.MessageMap[typing.Text, global___Symbol]: ...

    def __init__(self,
        *,
        symbols : typing.Optional[typing.Mapping[typing.Text, global___Symbol]] = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"symbols",b"symbols"]) -> None: ...
global___ListSymbolTableResult = ListSymbolTableResult

class Symbol(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    TYPE_FIELD_NUMBER: builtins.int
    type: typing.Text = ...

    def __init__(self,
        *,
        type : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"type",b"type"]) -> None: ...
global___Symbol = Symbol

class FunctionCallIdentifier(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    IDENTIFIER_FIELD_NUMBER: builtins.int
    identifier: typing.Text = ...

    def __init__(self,
        *,
        identifier : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"identifier",b"identifier"]) -> None: ...
global___FunctionCallIdentifier = FunctionCallIdentifier

class MessageFunctionArguments(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    NAME_FIELD_NUMBER: builtins.int
    name: typing.Text = ...

    def __init__(self,
        *,
        name : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"name",b"name"]) -> None: ...
global___MessageFunctionArguments = MessageFunctionArguments

class FunctionReturnValue(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    VALUE_FIELD_NUMBER: builtins.int
    value: typing.Text = ...

    def __init__(self,
        *,
        value : typing.Text = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"value",b"value"]) -> None: ...
global___FunctionReturnValue = FunctionReturnValue

class SetFunctionReturnValueResponse(google.protobuf.message.Message):
    DESCRIPTOR: google.protobuf.descriptor.Descriptor = ...
    RECEIVED_FIELD_NUMBER: builtins.int
    received: builtins.bool = ...

    def __init__(self,
        *,
        received : builtins.bool = ...,
        ) -> None: ...
    def ClearField(self, field_name: typing_extensions.Literal[u"received",b"received"]) -> None: ...
global___SetFunctionReturnValueResponse = SetFunctionReturnValueResponse

class SuperGSLCompiler(google.protobuf.service.Service, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def CreateCompilerSession(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___CreateCompilerSessionRequest,
        done: typing.Optional[typing.Callable[[global___CreateCompilerSessionResult], None]],
    ) -> concurrent.futures.Future[global___CreateCompilerSessionResult]: ...
    @abc.abstractmethod
    def DestroyCompilerSession(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___DestroyCompilerSessionRequest,
        done: typing.Optional[typing.Callable[[global___DestroyCompilerSessionResult], None]],
    ) -> concurrent.futures.Future[global___DestroyCompilerSessionResult]: ...
    @abc.abstractmethod
    def Compile(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___CompileRequest,
        done: typing.Optional[typing.Callable[[global___CompileResult], None]],
    ) -> concurrent.futures.Future[global___CompileResult]: ...
    @abc.abstractmethod
    def ListSymbolTable(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___ListSymbolTableRequest,
        done: typing.Optional[typing.Callable[[global___ListSymbolTableResult], None]],
    ) -> concurrent.futures.Future[global___ListSymbolTableResult]: ...
    @abc.abstractmethod
    def GetFunctionArguments(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___FunctionCallIdentifier,
        done: typing.Optional[typing.Callable[[global___MessageFunctionArguments], None]],
    ) -> concurrent.futures.Future[global___MessageFunctionArguments]: ...
    @abc.abstractmethod
    def SetFunctionReturnValue(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___FunctionReturnValue,
        done: typing.Optional[typing.Callable[[global___SetFunctionReturnValueResponse], None]],
    ) -> concurrent.futures.Future[global___SetFunctionReturnValueResponse]: ...
class SuperGSLCompiler_Stub(SuperGSLCompiler):
    def __init__(self, rpc_channel: google.protobuf.service.RpcChannel) -> None: ...
    def CreateCompilerSession(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___CreateCompilerSessionRequest,
        done: typing.Optional[typing.Callable[[global___CreateCompilerSessionResult], None]],
    ) -> concurrent.futures.Future[global___CreateCompilerSessionResult]: ...
    def DestroyCompilerSession(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___DestroyCompilerSessionRequest,
        done: typing.Optional[typing.Callable[[global___DestroyCompilerSessionResult], None]],
    ) -> concurrent.futures.Future[global___DestroyCompilerSessionResult]: ...
    def Compile(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___CompileRequest,
        done: typing.Optional[typing.Callable[[global___CompileResult], None]],
    ) -> concurrent.futures.Future[global___CompileResult]: ...
    def ListSymbolTable(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___ListSymbolTableRequest,
        done: typing.Optional[typing.Callable[[global___ListSymbolTableResult], None]],
    ) -> concurrent.futures.Future[global___ListSymbolTableResult]: ...
    def GetFunctionArguments(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___FunctionCallIdentifier,
        done: typing.Optional[typing.Callable[[global___MessageFunctionArguments], None]],
    ) -> concurrent.futures.Future[global___MessageFunctionArguments]: ...
    def SetFunctionReturnValue(self,
        rpc_controller: google.protobuf.service.RpcController,
        request: global___FunctionReturnValue,
        done: typing.Optional[typing.Callable[[global___SetFunctionReturnValueResponse], None]],
    ) -> concurrent.futures.Future[global___SetFunctionReturnValueResponse]: ...