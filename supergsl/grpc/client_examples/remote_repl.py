import grpc
import supergsl.grpc.stubs.sgsl_pb2_grpc as pb2_grpc
import supergsl.grpc.stubs.sgsl_pb2 as pb2


class SuperGSLCompilerClient(object):
    """Execute SuperGSL code through gRPC."""

    def __init__(self):
        self.host = 'localhost'
        self.server_port = 50051

        # instantiate a channel
        self.channel = grpc.insecure_channel(
            '{}:{}'.format(self.host, self.server_port))

        # bind the client and the server
        self.stub = pb2_grpc.SuperGSLCompilerStub(self.channel)

    def repl(self):
        """Await user input and then send it to a remote SuperGSL instance via gRPC.
        """
        result = self.stub.CreateCompilerSession(pb2.CreateCompilerSessionRequest())
        session_identifier = result.session_identifier

        while True:
            table_result = self.stub.ListSymbolTable(pb2.ListSymbolTableRequest(
                session_identifier=session_identifier
            ))
            print('Symbol Table')
            print(table_result)
            print('Sequence Store')
            sequence_store_result = self.stub.ListSequenceStore(
                pb2.ListSequenceStoreRequest(
                    session_identifier=session_identifier))

            print(sequence_store_result)

            user_input = input()

            request = pb2.CompileRequest(
                session_identifier=session_identifier,
                source_code=user_input)

            result = self.stub.Compile(request)
            print('RESULT', result.success, result.error)

if __name__ == '__main__':
    client = SuperGSLCompilerClient()
    client.repl()
