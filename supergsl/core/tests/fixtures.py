import mock
import random
from Bio.Seq import Seq
from supergsl.core.ast import Assembly, SymbolReference
from supergsl.types.builtin import PrimerPair
from supergsl.core.constants import THREE_PRIME
from supergsl.types.part import Part
from supergsl.types.position import SeqPosition
from supergsl.core.symbol_table import SymbolTable

class SuperGSLCoreFixtures(object):

    def mk_symbol_table(self):
        """Create a simple symbol table with a nested scope."""
        symbol_table = SymbolTable('awesome', None)
        symbol_table.insert('uHO', self.mk_part('uHO', 100)[1])

        nested_symbol_table = symbol_table.enter_nested_scope('awesome')
        nested_symbol_table.insert('tHUG', self.mk_part('tHUG', 33)[1])

        return symbol_table

    def get_assembly_ast(self):
        ast_part_nodes = [
            SymbolReference('uHO', None, False),
            SymbolReference('pADH1', None, False),
            SymbolReference('gERG10', None, False),
            SymbolReference('tADH1', None, False),
            SymbolReference('dHO', None, False),
        ]
        assembly = Assembly(ast_part_nodes)

        for part_node in ast_part_nodes:
            _, part_node.part = self.mk_part(part_node.identifier, 500, mk_primers=True)

        return assembly

    def mk_random_dna_sequence(self, seq_len):
        seq_str = ''.join(random.choice('CGTA') for _ in range(seq_len))
        return Seq(seq_str)

    def mk_extraction_primers(self, part):
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

    def mk_part(self, identifier, part_seq_len, mk_primers=True, roles=None):
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
