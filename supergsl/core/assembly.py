from typing import List
from supergsl.core.config import settings
from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.function import SuperGSLFunction
from supergsl.core.exception import ConfigurationException
from supergsl.utils import import_class


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

    def assemble(self):
        raise Exception('Not implemented. Subclass to implement.')


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assemblies):
        pass
