from typing import List
from Bio.Seq import Seq
from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.function import SuperGSLFunction
from supergsl.core.exception import ConfigurationException
from supergsl.utils import import_class


class Assembly(object):
    """Store the an assembled construct."""

    def get_sequence(self):
        """Return the complete sequence of the construct."""
        raise NotImplementedError('Subclass to implement.')

    def get_required_parts(self):
        """Return a list of parts required to construct this assembly."""
        raise NotImplementedError('Subclass to implement.')

    def get_part(self):
        """Retrieve a Part corresponding to this construct."""
        raise NotImplementedError('Subclass to implement.')


class AssemblerBase(SuperGSLFunction):

    # TODO(rmcl): How can we still retrieve some configuration for the assembler
    # from supergsl-config.json?

    #def __init__(self, name, config_options):
    #    self.name = name
    #    self.options = config_options

    def get_arguments(self):
        return []

    def get_return_type(self):
        return list

    def execute(self, sgsl_args, child_nodes=None):
        """
        """
        return self.assemble(child_nodes)

    def assemble(self, assemblies):
        raise Exception('Not implemented. Subclass to implement.')


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assemblies):
        """Iterate over each `ast.Assembly` node and generate an Assembly object."""
        for assembly_node in assemblies:
            assembly_sequence = Seq(''.join([
                part.get_sequence()
                for part in assembly_node.parts
            ]))

            assembly_node.assembly = Assembly(assembly_sequence)
