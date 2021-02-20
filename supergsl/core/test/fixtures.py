import mock
import random
from Bio.Seq import Seq
from supergsl.core.ast import Assembly, Part as AstPart
from supergsl.core.constants import THREE_PRIME
from supergsl.core.parts import Part, SeqPosition

class SuperGSLCoreFixtures(object):

    def get_assembly_ast(self):
        parts = [
            AstPart('uHO', None, False),
            AstPart('pADH1', None, False),
            AstPart('gERG10', None, False),
            AstPart('tADH1', None, False),
            AstPart('dHO', None, False),
        ]
        assembly = Assembly(parts)

        return assembly

    def mk_random_dna_sequence(self, seq_len):
        seq_str = ''.join(random.choice('CGTA') for _ in range(seq_len))
        return Seq(seq_str)

    def mk_part(self, identifier, part_seq_len):
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
            'TESTPART',
            start,
            end,
            provider=mock.Mock()
        )
        return reference_sequence, part
