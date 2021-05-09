from typing import List
from Bio.Seq import Seq

from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.function import SuperGSLFunction, SuperGSLFunctionDeclaration
from supergsl.core.types.part import Part
from supergsl.core.types.assembly import Assembly, AssemblyList


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
