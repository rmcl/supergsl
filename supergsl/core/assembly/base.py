from supergsl.core.config import settings
from supergsl.core.backend import BreadthFirstNodeFilteredPass
from supergsl.core.exception import ConfigurationException
from supergsl.utils import import_class


class AssemblerBase(object):

    def __init__(self, name, config_options):
        self.name = name
        self.options = config_options

    def assemble(self):
        raise Exception('Not implemented. Subclass to implement.')


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assemblies):
        pass

class AssemblerFactory(object):
    def __init__(self):
        self._initialize_assemblers()

    def get_assembler(self, assembler_name: str) -> AssemblerBase:
        try:
            return self._assemblers[assembler_name]
        except KeyError:
            raise Exception('Unknown assembler "%s".' % assembler_name)

    def _initialize_assemblers(self):
        self._assemblers : Dict[AssemblerBase] = {}

        if 'assemblers' not in settings:
            raise ConfigurationException('No assemblers have been defined. Check your supergGSL settings.')

        for assembler_config in settings['assemblers']:
            print('Initializing "%s"' % assembler_config['name'])
            assembler_class = import_class(assembler_config['assembler_class'])
            assembler_inst = assembler_class(
                assembler_config['name'],
                assembler_config['assembler_options'])

            self._assemblers[assembler_config['name']] = assembler_inst

class AssemblerPass(BreadthFirstNodeFilteredPass):
    """Visit each assembly block and execute the appropriate assembler."""

    def __init__(self):
        self._assembler_factory = AssemblerFactory()

    def get_node_handlers(self):
        return {
            'AssemblyBlock': self.visit_assembly_block_node,
        }

    def before_pass(self, ast):
        """Initialize the SBOL Document."""
        return ast

    def after_pass(self, ast):
        return ast

    def visit_assembly_block_node(self, node):
        assembler = self._assembler_factory.get_assembler(node.assembly_type)

        assembler.assemble(node.assemblies)
