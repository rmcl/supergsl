"""Implement the most basic assembler which concatenates a series of parts."""
from typing import List, Tuple
from supergsl.core.assembly import AssemblerBase, AssemblyLevel
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice
from supergsl.core.sequence import SliceMapping
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    Assembly,
    AssemblyResultSet
)


class FusionAssembler(AssemblerBase):
    """Create an assembly by fusing adjacent parts together without overlap.

    The most basic strategy is a "fusion" strategy where each part is annealed
    to its neighbors. The product of this strategy is likely useful for direct
    synthesis methods.

    .. code-block:: gsl

            from builtin import fuse
            from S288C import ADH1, TDH1, ERG10, HO

            fuse {
                HO_pADH1_gERG10: uHO ; pADH1 ; gERG10[1:728] ; dHO
                HO_pTDA1_gERG10: uHO ; pTDH1 ; gERG10[1:728] ; dHO
            }

    """

    def assemble_design(
        self,
        assembly_label : str,
        design_label : str,
        assembly_request : List[AssemblyLevel]
    ) -> Assembly:
        """Create an Assembly for a single design given a label and its parts."""
        design_label = f'ASM-{assembly_label}-{design_label}'

        cur_seq_pos = 0
        slice_mappings : List[SliceMapping] = []
        part_mappings : List[Tuple[Part, Slice]] = []
        for part in assembly_request:
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
        new_part = Part(
            design_label,
            assembly_sequence_entry,
            provider=self)

        assembly = Assembly(design_label, new_part, reagents=assembly_request)
        return assembly
