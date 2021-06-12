"""Implement the most basic assembler which concatenates a series of parts."""
from typing import List
from Bio.Seq import Seq
from supergsl.core.assembly import AssemblerBase
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    Assembly,
    AssemblyResultSet
)


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap."""

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyResultSet:
        """Iterate over `AssemblyDeclaration` and generate a set of assemblies."""

        assemblies : List[Assembly] = []
        for assembly_idx, assembly_request in enumerate(assembly_requests):

            designs = assembly_request.get_full_factorial_designs()
            for design_idx, design_parts in enumerate(designs):

                assembly_sequence = Seq(''.join([
                    str(part.sequence)
                    for part in design_parts
                ]))

                assembly_label = assembly_request.label or ('%03d' % assembly_idx)
                identifier = 'ASM-%s-%03d' % (
                    assembly_label,
                    design_idx)

                assembly = Assembly(identifier, assembly_sequence, design_parts)
                assemblies.append(assembly)

        return AssemblyResultSet(assemblies)
