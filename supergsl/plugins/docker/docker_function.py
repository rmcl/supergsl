import os
import inspect
from typing import Optional
import docker
from docker.errors import ImageNotFound
from supergsl.core.function import SuperGSLFunction
from supergsl.grpc import CompilerFunctionService

class DockerFunction(SuperGSLFunction):
    image_tag : Optional[str] = None

    def execute(self, params):
        """Invoke the function using docker."""

        print('Starting invoke of docker container function. %s' % self)
        grpc_server = self.start_grpc_server(
            self.serialize_params(params))

        command = self.get_docker_command(params)
        self.start_docker_container(command)

        # Wait for the functions to all return values or timeout.
        grpc_server.wait_for_function_completion(timeout=60)

        return grpc_server.get_return_value()

    def serialize_params(self, params : dict):
        return params

    def start_grpc_server(self, serialized_params : str):
        service = CompilerFunctionService(serialized_params)
        service.start_listening()
        return service

    def get_docker_command(self, params):
        """Generate the command string that will be passed to the docker container."""
        raise NotImplementedError('Subclass to implement.')

    def start_docker_container(self, command : str):
        """Run the docker container with the given command."""
        try:
            image = self.docker_client.images.get(self.get_image_tag())
        except ImageNotFound:
            image = self.build()

        container = self.docker_client.containers.run(
            image,
            command=command,
            auto_remove=True,
            detach=True)

        for log in container.logs(stream=True):
            print(log)


    @property
    def docker_client(self):
        return docker.from_env()


    def get_image_tag(self):
        if not self.image_tag:
            raise Exception('Image tag must be defined for docker plugin.')
        return self.image_tag


    def build(self):
        print('!--- Building Docker Image to invoke Plugin "%s" ----' % (
            self.__class__
        ))

        PLUGIN_BASE_DIR = os.path.dirname(
            os.path.abspath(inspect.getfile(self.__class__)))

        print(PLUGIN_BASE_DIR, 'tag', self.get_image_tag())

        build_log_gen = self.docker_client.api.build(
            path=PLUGIN_BASE_DIR,
            tag=self.get_image_tag(),
            pull=True,
            decode=True
        )

        self.print_docker_log(build_log_gen)

        print('!!--- Building Docker Image Complete "%s" ----' % (
            self.__class__
        ))

        return self.docker_client.images.get(self.get_image_tag())

    def print_docker_log(self, docker_log_gen):
        for chunk in docker_log_gen:
            if 'stream' in chunk:
                for line in chunk['stream'].splitlines():
                    print(line)
            elif 'status' in chunk:
                if chunk['status'] == 'Downloading':
                    print(chunk['progress'])
                else:
                    print(chunk, '????')
            else:
                print(chunk, 'NOT A STREAM')
