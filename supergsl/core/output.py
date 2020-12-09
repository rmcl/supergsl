import csv
from typing import Optional
from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import ConfigurationException
from supergsl.core.config import settings
from supergsl.utils import import_class


class OutputProvider(BreadthFirstNodeFilteredPass):
    name : Optional[str] = None

    @classmethod
    def get_output_name(cls):
        if cls.name is None:
            raise ConfigurationException('%s does not specify its output name.')

        return cls.name

class ASTPrintOutputProvider(OutputProvider):
    name = 'print'

    def before_pass(self, ast):
        """Initialize the SBOL Document."""

        import pprint
        pprint.pprint(ast.eval())

        return ast

class PrimerOutputProvider(OutputProvider):
    name = 'primers'

    def get_node_handlers(self):
        return {
            'Part': self.visit_part_node,
        }

    def visit_part_node(self, node):
        part = node.part
        self.primers[part.identifier] = {
            'Forward Primer': part.forward_primer,
            'Reverse Primer': part.reverse_primer
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        self.primers = {}
        return ast

    def after_pass(self, ast):
        """Save a TSV of part primers."""

        with open('primers.txt', 'w+') as fp:
            writer = csv.DictWriter(fp, (
                'Part Identifier',
                'Forward Primer',
                'Reverse Primer'
            ))

            writer.writeheader()

            for part_identifier, details in self.primers.items():
                output = details.copy()
                output['Part Identifier'] = part_identifier

                writer.writerow(output)

        return ast


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
        return ast

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

    def __init__(self):
        self.resolve_providers()

    def validate_args(self, args):
        self.desired_output_providers = []
        for output_format_name in args.output_format:
            try:
                outputer_class = self.available_outputers[output_format_name]
            except KeyError:
                raise Exception('Unknown output format "%s".' % output_format_name)

            output_inst = outputer_class(None)
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
            print('running %s' % output_pass_inst.name)
            output_pass_inst.perform(ast)
