import mock
import random
from Bio.Seq import Seq
from supergsl.core.ast import Assembly, Part as AstPart
from supergsl.core.types import PrimerPair
from supergsl.core.constants import THREE_PRIME
from supergsl.core.parts import Part, SeqPosition

class SuperGSLCoreFixtures(object):

    def get_assembly_ast(self):
        ast_part_nodes = [
            AstPart('uHO', None, False),
            AstPart('pADH1', None, False),
            AstPart('gERG10', None, False),
            AstPart('tADH1', None, False),
            AstPart('dHO', None, False),
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

    def mk_part(self, identifier, part_seq_len, mk_primers=True):
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
            provider=mock.Mock()
        )

        if mk_primers:
            primer_pair = self.mk_extraction_primers(part)
            part.set_extraction_primers(primer_pair)

        return reference_sequence, part
