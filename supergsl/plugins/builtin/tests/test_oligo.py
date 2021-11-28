"""Unit tests for the Fuse Assembler."""
import unittest
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.types.assembly import AssemblyDeclaration
from supergsl.plugins.builtin.oligos import SyntheticOligoAssembler


class SyntheticOligoAssemblerTestCase(unittest.TestCase):
    """Test that the parser correctly tokens into a valid AST."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.oligo_overlap = 20
        self.assembler = SyntheticOligoAssembler(self.fixtures.mk_function_config_object({
            'max_oligo_len': 45,
            'max_num_oligos': 5,
            'min_overlap_len': self.oligo_overlap
        }))

    def test_assemble_concatenates_part_sequences(self):
        """Fuse assembler is expected to append part sequences together."""
        five_prime_flank = self.fixtures.mk_part('5p_flank', 25)[1]
        promoter = self.fixtures.mk_part('promoter', 30)[1]
        three_prime_flank = self.fixtures.mk_part('3p_flank', 25)[1]

        declaration = AssemblyDeclaration('SweetAssembly', [
            five_prime_flank,
            promoter,
            three_prime_flank
        ])

        assembly_result_set = list(self.assembler.assemble([declaration]))
        assemblies = list(assembly_result_set)
        self.assertEqual(len(assemblies), 1)
        self.assertEqual(assemblies[0].sequence, Seq(''.join([
            str(five_prime_flank.sequence),
            str(promoter.sequence),
            str(three_prime_flank.sequence)
        ])))

        primers = assemblies[0].reagents

        expected_primer_sequences = []
        for primer_idx, primer in enumerate(primers):
            if primer_idx == 0:
                expected_primer_sequences.append(str(primer.sequence))
            else:
                expected_primer_sequences.append(str(primer.sequence[self.oligo_overlap:]))

        expected_part_sequence = ''.join(expected_primer_sequences)
        self.assertEqual(expected_part_sequence, assemblies[0].sequence)

    '''
    def test_assemble_concatenates_part_sequences_with_collection(self):
        """Fuse assembler is expected to append part sequences together."""
        five_prime_flank = self.fixtures.mk_part('5p_flank', 25)[1]
        promoter_collection = self.fixtures.mk_part_collection(num_parts=3, part_len=20)
        three_prime_flank = self.fixtures.mk_part('3p_flank', 25)[1]

        declaration = AssemblyDeclaration('SweetAssembly', [
            five_prime_flank,
            promoter_collection,
            three_prime_flank
        ])


        result = list(self.assembler.assemble([declaration]))

        assemblies = list(result)
        self.assertEqual(assemblies[0].sequence, 'HI')
    '''
