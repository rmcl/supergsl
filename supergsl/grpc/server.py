from typing import Optional
import grpc
from grpc import Server
from concurrent import futures
import supergsl.grpc.stubs.sgsl_pb2_grpc as pb2_grpc
import supergsl.grpc.stubs.sgsl_pb2 as pb2


class CompilerFunctionService(pb2_grpc.CompilerFunctionServicer):

    def __init__(self, function_invoke_params : dict):
        self.function_invoke_params : dict = function_invoke_params
        self.server : Optional[Server] = None

    def SetFunctionReturnValue(self, request, context):
        """Missing associated documentation comment in .proto file."""

        message = request.value
        print('message')
        return pb2.SetFunctionReturnValueResponse(received=True)

    def start_listening(self):
        """Start the gRPC server in a new thread."""
        self.server = grpc.server(futures.ThreadPoolExecutor(max_workers=1))
        pb2_grpc.add_CompilerFunctionServicer_to_server(self, self.server)
        self.server.add_insecure_port('[::]:50051')
        self.server.start()

    def wait_for_function_completion(self, timeout = None):
        """Wait for the invoked functions to set their return value.

        If timeout ellapsed then terminate server and raise an exception.
        """
        self.server.wait_for_termination()
