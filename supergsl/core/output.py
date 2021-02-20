import csv
import pprint
from typing import Optional
from supergsl.core.ast import Node
from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import ConfigurationException
from supergsl.utils import import_class


class OutputProvider(BreadthFirstNodeFilteredPass):
    """OutputProviders perform an AST pass to generate a usable output.
    This is a base class."""

    name: Optional[str] = None

    @classmethod
    def get_output_name(cls):
        """Get the Output Provider's name."""
        if cls.name is None:
            raise ConfigurationException('%s does not specify its output name.')

        return cls.name

    def before_pass(self, ast : Node) -> Node:
        pass

    def after_pass(self, ast : Node) -> Node:
        pass

class ASTPrintOutputProvider(OutputProvider):
    """Generate a pretty print output of the AST and write it to stdout."""
    name = 'print'
    stream = None

    def before_pass(self, ast):
        """Initialize the SBOL Document."""

        printer = pprint.PrettyPrinter(stream=self.stream)
        printer.pprint(ast.eval())


class PrimerOutputProvider(OutputProvider):
    name = 'primers'

    def get_node_handlers(self):
        return {
            'Part': self.visit_part_node,
        }

    def visit_part_node(self, node):
        """Visit each part node and record the part's extraction primers."""
        part = node.part
        self.primers[part.identifier] = {
            'Part Identifier': part.identifier,
            'Forward Primer': str(part.extraction_primers.forward_primer.get_sequence()),
            'Reverse Primer': str(part.extraction_primers.reverse_primer.get_sequence())
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        self.primers = {}

    def _open_primer_file(self):
        return open('primers.txt', 'w+')

    def after_pass(self, ast):
        """Save a TSV of part primers."""

        print(self.primers)
        with self._open_primer_file() as fp:
            writer = csv.DictWriter(fp, (
                'Part Identifier',
                'Forward Primer',
                'Reverse Primer'
            ))

            writer.writeheader()

            for _, details in self.primers.items():
                output = details.copy()
                writer.writerow(output)



class TestOutputProvider(OutputProvider):
    """Output provider used by Integration tests to introspect the output of the compiler."""
    name = 'test'

    def get_node_handlers(self):
        return {
            'Assembly': self.visit_assembly_node,
            'Part': self.visit_part_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        self.parts = []
        self.assemblies = []

    def visit_part_node(self, node):
        self.parts.append(node.part)

    def visit_assembly_node(self, node):
        assembly_parts = [
            part.identifier
            for part in node.parts
        ]

        assembly_idx = len(self.assemblies)
        self.assemblies.append({
            'identifier': assembly_idx,
            'parts': assembly_parts
        })

    def get_assemblies(self):
        return self.assemblies

    def get_parts(self):
        return self.parts



class OutputPipeline(object):
    """Store the results of the compilation in user specified output formats."""

    def __init__(self, compiler_settings):
        self.settings = compiler_settings
        self.resolve_providers()

    def validate_args(self, args):
        self.desired_output_providers = []
        for output_format_name in args.output_format:
            try:
                outputer_class = self.available_outputers[output_format_name]
            except KeyError:
                raise Exception('Unknown output format "%s".' % output_format_name)

            output_inst = outputer_class(None, allow_modification=False)
            self.desired_output_providers.append(output_inst)

    def resolve_providers(self):
        """Resolve the output providers from the supergsl-config settings."""

        self.available_outputers = {}
        for provider_class_path in self.settings['output_providers']:
            provider_class = import_class(provider_class_path)

            if not issubclass(provider_class, OutputProvider):
                raise ConfigurationException(
                    '"%s" is not an instance of OutputProvider' % provider_class_path)

            self.available_outputers[provider_class.get_output_name()] = provider_class


    def run(self, ast, args):
        """Output the compiled AST into the desired output formats."""

        for output_pass_inst in self.desired_output_providers:
            print('running %s' % output_pass_inst.name)
            output_pass_inst.perform(ast)
