"""Define fixture data for testing SuperGSL."""
from unittest.mock import Mock
import random
from typing import Tuple, List
from Bio.Seq import Seq
from supergsl.core.types.primer import PrimerPair
from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.part import Part
from supergsl.core.types.builtin import Collection
from supergsl.core.types.position import SeqPosition
from supergsl.core.types.assembly import (
    Assembly,
    AssemblyDeclaration,
    AssemblyResultSet
)
from supergsl.core.symbol_table import SymbolTable

class SuperGSLCoreFixtures(object):

    def mk_symbol_table(self) -> SymbolTable:
        """Create a simple symbol table with a nested scope."""
        symbol_table = SymbolTable('awesome', None)
        symbol_table.insert('uHO', self.mk_part('uHO', 100)[1])

        nested_symbol_table = symbol_table.enter_nested_scope('awesome')
        nested_symbol_table.insert('tHUG', self.mk_part('tHUG', 33)[1])

        return symbol_table

    def mk_random_dna_sequence(self, seq_len) -> Seq:
        seq_str = ''.join(random.choice('CGTA') for _ in range(seq_len))
        return Seq(seq_str)

    def mk_extraction_primers(self, part) -> PrimerPair:
        """Primers are often complementary sequences flanking the DNA sequence of interest.

        Real primers woud use some smarter algorithm to make sure the resulting
        DNA sequence has a melting temperature that is conducive to the parameters
        of a the particular PCR thermocycler.
        """
        complement_sequence = part.sequence.complement()
        return PrimerPair.from_sequences(
            complement_sequence[:20],
            complement_sequence[-20:]
        )

    def mk_part(
        self,
        identifier,
        part_seq_len,
        mk_primers=True,
        roles=None
    ) -> Tuple[Seq,Part]:
        """Create a mock Part.

        Part is derived from a reference sequence three times longer than the part
        itself.
        """
        ref_seq_len = part_seq_len * 3
        reference_sequence = self.mk_random_dna_sequence(ref_seq_len)

        start = SeqPosition.from_reference(
            x=part_seq_len,
            rel_to=THREE_PRIME,
            approximate=False,
            reference=reference_sequence
        )
        end = start.get_relative_position(x=part_seq_len)

        part = Part(
            identifier,
            start,
            end,
            provider=Mock(),
            roles=roles
        )

        if mk_primers:
            primer_pair = self.mk_extraction_primers(part)
            part.set_extraction_primers(primer_pair)

        return reference_sequence, part

    def mk_part_collection(self, num_parts=3) -> Collection:
        """Return a `Collection` of parts."""

        return Collection([
            self.mk_part('pGAL%d' % index, 100)[1]
            for index in range(num_parts)
        ])

    def mk_assembly_declaration_ex1(self, name='test_declaration') -> AssemblyDeclaration:
        """Make a AssemblyDeclaration with 4 parts including one with a collection of three parts"""
        promoter_collection = self.mk_part_collection(num_parts=3)
        gene = self.mk_part('gGENE', 500)[1]
        upstream = self.mk_part('uHO', 50)[1]
        downstream = self.mk_part('dHO', 50)[1]

        return AssemblyDeclaration(name, [
            upstream,
            promoter_collection,
            gene,
            downstream
        ])

    def mk_assembly(self, identifier='asm1', num_parts=2) -> Assembly:
        """Create a `Assembly` containing num_parts with random sequences of len 0 to 100."""
        parts : List[Part] = list([
            self.mk_part('part-%03d' % part_index, random.randint(20, 100))[1]
            for part_index in range(num_parts)
        ])
        sequence = Seq(''.join([
            str(part.sequence)
            for part in parts
        ]))
        return Assembly(identifier, sequence, parts)

    def mk_assembly_result_set(self, num_assembly=2) -> AssemblyResultSet:
        """Create an AssemblyResultSet with `num_assembly` assemblies."""
        return AssemblyResultSet([
            self.mk_assembly(
                'asm-%d' % assembly_index,
                num_parts=random.randint(1,5))
            for assembly_index in range(num_assembly)
        ])
