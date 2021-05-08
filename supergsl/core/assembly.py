from typing import List
from Bio.Seq import Seq

from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.exception import ConfigurationError
from supergsl.core.types import SuperGSLType
from supergsl.core.parts import Part
from supergsl.utils import import_class


class Assembly(SuperGSLType):
    """Store the an assembled construct."""

    def __init__(self, identifier : str, sequence : Seq, parts : List[Part]):
        self.identifier = identifier
        self.sequence = sequence
        self.parts = parts

    def get_identifier(self):
        return self.identifier

    def get_sequence(self):
        """Return the complete sequence of the construct."""
        return self.sequence

    def get_required_parts(self):
        """Return a list of parts required to construct this assembly."""
        return self.parts

    def get_part(self):
        """Retrieve a Part corresponding to this construct."""
        raise NotImplementedError('Subclass to implement.')


class AssemblyList(SuperGSLType):
    """A collection of Assemblies."""
    def __init__(self, assemblies : Assembly):
        self.assemblies = assemblies

    def __iter__(self):
        return iter(self.assemblies)

class AssemblerBase(SuperGSLFunction):
    """Base class for functions implementing Assemblers."""
    def get_arguments(self):
        return []

    def get_return_type(self):
        return AssemblyList

    def execute(self, params):
        return self.assemble(params['children'])

    def assemble(self, assembly_part_lists : List[List[Part]]) -> AssemblyList:
        """Iterate over `Part` and generate an Assembly object."""
        raise NotImplementedError('Not implemented. Subclass to implement.')


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assembly_part_lists : List[List[Part]]) -> AssemblyList:
        """Iterate over `Part` and generate an Assembly object."""

        assemblies = []
        for assembly_idx, assembly_part_list in enumerate(assembly_part_lists):
            assembly_sequence = Seq(''.join([
                str(part.get_sequence().seq)
                for part in assembly_part_list
            ]))

            identifier = str('ASM-%05d' % assembly_idx)
            assembly = Assembly(identifier, assembly_sequence, assembly_part_list)
            assemblies.append(assembly)

        return AssemblyList(assemblies)


class BuiltinAssemblersPlugin(SuperGSLPlugin):
    """Plugin stub to help register basic Assemblers."""

    def register(self, symbol_table, compiler_settings):
        """Register built in assemblers."""
        symbol_table.insert('fuse',
            SuperGSLFunctionDeclaration(FusionAssembler, compiler_settings))
