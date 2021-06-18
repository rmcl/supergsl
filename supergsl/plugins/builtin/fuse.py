"""Implement the most basic assembler which concatenates a series of parts."""
from typing import List
from Bio.Seq import Seq
from supergsl.core.constants import THREE_PRIME
from supergsl.core.assembly import AssemblerBase
from supergsl.core.types.part import Part
from supergsl.core.types.position import SeqPosition
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

                assembly = Assembly(identifier, assembly_sequence)
                cur_seq_pos = 0
                for part in design_parts:
                    start_pos, end_pos = self.get_part_position(
                        assembly_sequence,
                        part,
                        cur_seq_pos)

                    assembly.add_part(part, start_pos, end_pos)
                    cur_seq_pos += len(part.sequence)
                assemblies.append(assembly)

        return AssemblyResultSet(assemblies)

    def get_part_position(self, assembly_sequence : Seq, part : Part, start_pos : int):
        """Create SeqPosition objects in the provided assembly sequence."""
        start = SeqPosition.from_reference(
            x=start_pos,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=assembly_sequence
        )
        end = start.get_relative_position(len(part.sequence))
        return start, end
