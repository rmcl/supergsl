import os
import inspect
import docker
from docker.errors import ImageNotFound
from supergsl.core.function import SuperGSLFunction


class DockerFunction(SuperGSLFunction):
    image_tag = None

    def execute(self, sgsl_args, child_nodes=None):
        """Invoke the function using docker."""

        # Setup gRPC????
        # Invoke the docker container

        print('CUT IT UP!')


        self.invoke()

        return None



    @property
    def docker_client(self):
        return docker.from_env()

    def get_image_tag(self):
        if not self.image_tag:
            raise Exception('Image tag must be defined for docker plugin.')
        return self.image_tag

    def invoke(self):
        try:
            image = self.docker_client.images.get(self.get_image_tag())
        except ImageNotFound:
            image = self.build()

        container = self.docker_client.containers.run(
            image,
            command='hello',
            auto_remove=True,
            detach=True)

        for log in container.logs(stream=True):
            print(log)

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
