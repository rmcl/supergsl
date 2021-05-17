from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunctionDeclaration
from supergsl.core.types.builtin import NucleotideSequence
from supergsl.plugins.docker import DockerFunction

class ExampleDockerFunction(DockerFunction):
    """An example of a Docker wrapped SuperGSLFunction.

    Return a random nucleotide sequence generated within a
    python script executed in a docker container.
    """

    image_tag = 'example_docker:latest'
    name = 'random_sequence'

    def get_arguments(self):
        return [
            ('num_results', int)
        ]

    def get_return_type(self):
        return NucleotideSequence


class ExampleDockerPlugin(SuperGSLPlugin):
    """Register example docker plugin."""

    def register(self, compiler_settings):
        """Register functions provide by chopchop."""
        self.symbol_table.register(
            'examples',
            'docker_example',
            SuperGSLFunctionDeclaration(
                ExampleDockerFunction,
                compiler_settings))
