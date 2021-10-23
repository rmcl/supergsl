"""Implement the most basic assembler which concatenates a series of parts."""
from typing import List, Tuple
from Bio.Seq import Seq
from supergsl.core.constants import THREE_PRIME
from supergsl.core.assembly import AssemblerBase
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice
from supergsl.core.sequence import SliceMapping
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

                assembly_label = assembly_request.label or ('%03d' % assembly_idx)
                design_label = 'ASM-%s-%03d' % (
                    assembly_label,
                    design_idx)

                assemblies.append(
                    self.assemble_one_design(design_label, design_parts))

        return AssemblyResultSet(assemblies)

    def assemble_one_design(self, design_label, design_parts):
        """Create an Assembly for a single design given a label and its parts."""
        cur_seq_pos = 0
        slice_mappings : List[SliceMapping] = []
        part_mappings : List[Tuple[Part, Slice]] = []
        for part in design_parts:
            start_pos = cur_seq_pos
            end_pos = cur_seq_pos + len(part.sequence)
            cur_seq_pos = end_pos

            source_slice = Slice.from_entire_sequence()
            target_slice = Slice.from_five_prime_indexes(
                start_index=start_pos,
                end_index=end_pos)

            slice_mappings.append(
                SliceMapping(part.sequence_entry, source_slice, target_slice))
            part_mappings.append(
                (part, target_slice))

        assembly_sequence_entry = self.sequence_store.concatenate(slice_mappings)

        assembly = Assembly(design_label, assembly_sequence_entry)
        for part, target_slice in part_mappings:
            assembly.add_part(part, target_slice)

        return assembly
