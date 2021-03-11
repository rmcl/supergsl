from supergsl.core.plugin import SuperGSLPlugin
from supergsl.plugins.docker import DockerFunction

class ExampleDockerFunction(DockerFunction):
    """An example of a Docker wrapped SuperGSLFunction.

    Return a random nucleotide sequence generated within a
    python script executed in a docker container.
    """

    name = 'random_sequence'

    def get_arguments(self):
        return [
            argument('num_results', int)
        ]

    def get_return_type(self):
        return NucleotideConstant


class ExampleDockerPlugin(SuperGSLPlugin):

    def register(self, symbol_table, compiler_settings):
        """Register functions provide by chopchop."""
        symbol_table.register('examples', ExampleDockerPlugin())
