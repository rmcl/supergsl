"""Implement a Assembler useful for creating parts via overlapping short oligonucleotides."""
from typing import List
from Bio.Seq import Seq
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    Assembly,
    AssemblyResultSet
)
from supergsl.core.assembly import AssemblerBase

DEFAULT_MAX_OLIGO_LEN = 200
DEFAULT_MIN_OVERLAP = 20
DEFAULT_NUM_OLIGOS = 2


class SyntheticOligoAssembler(AssemblerBase):
    """Create parts to be constructed by ordering a set of short synthetic oligos.

    These oligos are typically used as primers, but can be used for short parts
    without a template.
    """

    def __init__(self, config_options):
        self.max_oligo_len = config_options.get('max_oligo_len', DEFAULT_MAX_OLIGO_LEN)
        self.min_overlap = config_options.get('max_oligo_len', DEFAULT_MIN_OVERLAP)
        self.max_num_oligos = config_options.get('max_num_oligos', DEFAULT_NUM_OLIGOS)

    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyResultSet:
        """Iterate over `Part` and generate an Assembly object."""

        oligos : List[Seq] = []

        assemblies : List[Assembly] = []
        for assembly_idx, assembly_request in enumerate(assembly_requests):
            parts = assembly_request.get_levels_by_factor_type('Part')

            assembly_sequence = Seq(''.join([
                str(part.get_sequence().seq)
                for part in parts
            ]))

            oligo_idx = 0
            seq_pos = 0
            sequence_length = len(assembly_sequence)
            while True:
                if oligo_idx == self.max_num_oligos:
                    raise Exception((
                        'Assembly is to large to be synthesized. Construct length: {} '
                        'Max Oligo Length: {} Max Oligos: {} ').format(
                            len(assembly_sequence),
                            self.max_oligo_len,
                            self.max_num_oligos
                        ))

                remaining_seq_len = sequence_length - seq_pos
                if remaining_seq_len < self.max_oligo_len:
                    seq_pos = min(0, seq_pos - (self.max_oligo_len - remaining_seq_len))

                oligos.append(assembly_sequence[seq_pos:self.max_oligo_len])
                seq_pos += self.max_oligo_len
                if seq_pos >= sequence_length:
                    break

                seq_pos -= self.min_overlap
                oligo_idx += 1

            print(oligos)

            identifier = str('ASM-%05d' % assembly_idx)
            assembly = Assembly(identifier, assembly_sequence, oligos)
            assemblies.append(assembly)

        return AssemblyResultSet(assemblies)
