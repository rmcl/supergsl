"""Implement a Assembler useful for creating parts via overlapping short oligonucleotides."""
from typing import List
from Bio.Seq import Seq
from supergsl.core.function import SuperGSLFunctionConfig
from supergsl.core.sequence import SequenceEntry, SliceMapping, Role
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice
from supergsl.core.types.primer import Primer
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    Assembly,
    AssemblyResultSet
)
from supergsl.core.assembly import AssemblerBase

DEFAULT_MAX_OLIGO_LEN = 200
DEFAULT_MIN_OVERLAP = 20
DEFAULT_NUM_OLIGOS = 2


synthetic_oligo_role = Role(
    uri='builtin/synthetic-oligo',
    name='Synthetic Oligo',
    description='A short synthesized piece of DNA.'
)


class SyntheticOligoAssembler(AssemblerBase):
    """Create parts to be constructed by ordering a set of short synthetic oligos.

    These oligos are typically used as primers, but can be used for short parts
    without a template.

    Example usage

    In the iGEM community, there are some small parts such as bacterial promoters
    where it is more convenient (and cheaper) to simply order a set of short
    synthetic oligos. In this example we prepend and append the biobrick prefix
    sequences to a bacterial promoter from `Anderson`. The resulting part can
    then be used using the biobick 3a assembly method to make more interesting
    and larger parts.

    Given this snippet of SuperGSL code:

    .. code-block:: gsl

        # http://parts.igem.org/Help:Promoters/Construction
        # http://parts.igem.org/Promoters/Catalog/Anderson

        from biobrick import bbPrefix, bbSuffix

        # http://parts.igem.org/Part:BBa_J23100
        from synbiohub import BBa_J23100

        synthetic-oligos {
            bbPrefix ; BBa_J23100 ; bbSuffix
        }

    As a further example, if you wanted to create a whole collection of biobrick
    promoter parts you could use SuperGSL collection syntax to generate the full
    factorial of promoters with flanking biobrick regions.


    .. code-block:: gsl

        from biobrick import bbPrefix, bbSuffix
        from anderson_promoters import promoters

        synthetic-oligos {
            bbPrefix ; promoters ; bbSuffix
        }

    """

    def __init__(self, config : SuperGSLFunctionConfig):
        config_options = config.compiler_settings
        self.sequence_store = config.sequence_store
        self.max_oligo_len = config_options.get('max_oligo_len', DEFAULT_MAX_OLIGO_LEN)
        self.min_overlap = config_options.get('min_overlap_len', DEFAULT_MIN_OVERLAP)
        self.max_num_oligos = config_options.get('max_num_oligos', DEFAULT_NUM_OLIGOS)


    def build_one_sequence_entry(self, design_parts : List[Part]):
        """Create a SequenceEntry for a single design given its parts."""
        cur_seq_pos = 0
        slice_mappings : List[SliceMapping] = []

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

        return self.sequence_store.concatenate(slice_mappings)

    def build_oligos_for_sequence_entry(self, sequence_entry : SequenceEntry) -> List[Primer]:
        sequence_length = sequence_entry.sequence_length

        oligo_idx = 0
        cur_seq_pos = 0
        oligo_entries : List[Primer] = []
        while True:
            if oligo_idx == self.max_num_oligos:
                raise Exception((
                    'Assembly is to large to be synthesized. Construct length: {} '
                    'Max Oligo Length: {} Max Oligos: {} ').format(
                        sequence_length,
                        self.max_oligo_len,
                        self.max_num_oligos
                    ))

            remaining_seq_len = sequence_length - cur_seq_pos
            if remaining_seq_len < self.max_oligo_len:
                end_pos = cur_seq_pos + remaining_seq_len
            else:
                end_pos = cur_seq_pos + self.max_oligo_len

            oligo_entry = self.sequence_store.slice(
                sequence_entry,
                Slice.from_five_prime_indexes(cur_seq_pos, end_pos),
                new_sequence_roles = [synthetic_oligo_role],
                entry_link_roles = [synthetic_oligo_role]) # TODO may need a different role for the entry link.

            new_primer = Primer(oligo_entry)
            oligo_entries.append(new_primer)

            cur_seq_pos = end_pos
            if cur_seq_pos >= sequence_length:
                break

            cur_seq_pos -= self.min_overlap
            oligo_idx += 1

        return oligo_entries


    def assemble(self, assembly_requests : List[AssemblyDeclaration]) -> AssemblyResultSet:
        """Iterate over `Part` and generate an Assembly object."""

        oligos : List[Seq] = []
        assemblies : List[Assembly] = []

        for assembly_idx, assembly_request in enumerate(assembly_requests):
            designs = assembly_request.get_full_factorial_designs()
            for design_idx, design_parts in enumerate(designs):
                sequence_entry = self.build_one_sequence_entry(design_parts)
                oligos = self.build_oligos_for_sequence_entry(sequence_entry)

                identifier = str('ASM-%03d-%05d' % (design_idx, assembly_idx))

                new_part = Part(
                    identifier,
                    sequence_entry,
                    provider=self,
                )

                assembly = Assembly(
                    identifier,
                    part=new_part,
                    reagents=oligos)
                assemblies.append(assembly)

        return AssemblyResultSet(assemblies)
