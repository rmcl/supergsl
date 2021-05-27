import grpc
import supergsl.grpc.stubs.sgsl_pb2_grpc as pb2_grpc
import supergsl.grpc.stubs.sgsl_pb2 as pb2


class SuperGSLFunctionInvocationClient(object):
    """Allow a remote SuperGSLFunction to access details of its invocation.

    This is useful for Docker based plugins so they can interact with the compiler
    which is executing on the host machine or in a different container.
    """

    def __init__(self):
        self.host = 'localhost'
        self.server_port = 50051

        # instantiate a channel
        self.channel = grpc.insecure_channel(
            '{}:{}'.format(self.host, self.server_port))

        # bind the client and the server
        self.stub = pb2_grpc.SuperGSLCompilerStub(self.channel)

    def get_params_and_update_return_val(self):
        """
        Client function to call the rpc for GetServerResponse
        """


if __name__ == '__main__':
    client = SuperGSLFunctionInvocationClient()
    client.get_params_and_update_return_val()
