"""Implement the most basic assembler which concatenates a series of parts."""
from typing import List
from Bio.Seq import Seq
from supergsl.core.assembly import AssemblerBase
from supergsl.types.assembly import (
    AssemblyDeclaration,
    Assembly,
    AssemblyList
)


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyList:
        """Iterate over `Part` and generate an Assembly object."""

        assemblies : List[Assembly] = []
        for assembly_idx, assembly_request in enumerate(assembly_requests):
            parts = assembly_request.get_parts()
            assembly_sequence = Seq(''.join([
                str(part.get_sequence().seq)
                for part in parts
            ]))

            identifier = str('ASM-%05d' % assembly_idx)
            assembly = Assembly(identifier, assembly_sequence, parts)
            assemblies.append(assembly)

        return AssemblyList(assemblies)
