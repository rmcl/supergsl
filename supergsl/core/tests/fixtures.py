import mock
import random
from typing import Tuple
from Bio.Seq import Seq

from supergsl.core.constants import THREE_PRIME
from supergsl.core.types.assembly import AssemblyDeclaration
from supergsl.core.symbol_table import SymbolTable

from supergsl.types.builtin import PrimerPair, Collection
from supergsl.types.part import Part
from supergsl.types.position import SeqPosition


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
        complement_sequence = part.get_sequence().seq.complement()
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
            provider=mock.Mock(),
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
