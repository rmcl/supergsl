from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import ConfigurationException
from supergsl.core.config import settings
from supergsl.utils import import_class


class OutputProvider(BreadthFirstNodeFilteredPass):
    name = None

    @classmethod
    def get_output_name(cls):
        if cls.name is None:
            raise ConfigurationException('%s does not specify its output name.')

        return cls.name


class OutputPipeline(object):
    """Store the results of the compilation in user specified output formats."""

    def __init__(self):
        self.resolve_providers()

    def validate_args(self, args):
        self.desired_output_providers = []
        for output_format_name in args.output_format:
            try:
                outputer_class = self.available_outputers[output_format_name]
            except KeyError:
                raise Exception('Unknown output format "%s".' % output_format_name)

            output_inst = outputer_class()
            self.desired_output_providers.append(output_inst)

    def resolve_providers(self):
        """Resolve the output providers from the supergsl-config settings."""

        self.available_outputers = {}
        for provider_class_path in settings['output_providers']:
            provider_class = import_class(provider_class_path)

            if not issubclass(provider_class, OutputProvider):
                raise ConfigurationException(
                    '"%s" is not an instance of OutputProvider' % provider_class_path)

            self.available_outputers[provider_class.get_output_name()] = provider_class


    def run(self, ast, args):
        """Output the compiled AST into the desired output formats."""

        for output_pass_inst in self.desired_output_providers:
            output_pass_inst.perform(ast)
