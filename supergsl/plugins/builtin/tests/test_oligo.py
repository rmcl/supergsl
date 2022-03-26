"""Unit tests for the Fuse Assembler."""
import unittest
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.types.assembly import (
    AssemblyDeclaration,
    AssemblyLevelDeclaration
)
from supergsl.plugins.builtin.oligos import SyntheticOligoAssembler


class SyntheticOligoAssemblerTestCase(unittest.TestCase):
    """Test that the parser correctly tokens into a valid AST."""
    maxDiff = None

    def setUp(self):
        self.fixtures = SuperGSLCoreFixtures()
        self.oligo_overlap = 20
        self.assembler = SyntheticOligoAssembler(self.fixtures.mk_provider_config({
            'max_oligo_len': 80,
            'max_num_oligos': 5,
            'min_overlap_len': self.oligo_overlap
        }))

    def test_assemble_concatenates_part_sequences(self):
        """Fuse assembler is expected to append part sequences together."""
        five_prime_flank = self.fixtures.mk_part('5p_flank', 25)[1]
        promoter = self.fixtures.mk_part('promoter', 30)[1]
        three_prime_flank = self.fixtures.mk_part('3p_flank', 25)[1]

        declaration = AssemblyDeclaration('SweetAssembly', [
            AssemblyLevelDeclaration(five_prime_flank, None),
            AssemblyLevelDeclaration(promoter, None),
            AssemblyLevelDeclaration(three_prime_flank, None)
        ])

        assembly_result_set = list(self.assembler.execute({
            'children': [declaration],
            'constraints': []
        }))
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


    def test_assemble_concatenates_part_sequences_with_collection(self):
        """Fuse assembler is expected to append part sequences together also if their is a collection."""
        five_prime_flank = self.fixtures.mk_part('5p_flank', 25)[1]
        promoter_collection = self.fixtures.mk_part_collection(num_parts=3, part_len=20)
        three_prime_flank = self.fixtures.mk_part('3p_flank', 25)[1]

        declaration = AssemblyDeclaration('SweetAssembly', [
            AssemblyLevelDeclaration(five_prime_flank, None),
            AssemblyLevelDeclaration(promoter_collection, None),
            AssemblyLevelDeclaration(three_prime_flank, None)
        ])

        result = list(self.assembler.execute({
            'children': [declaration],
            'constraints': []
        }))

        assemblies = list(result)

        # 1 upstream; 3 parts in collection; 1 downstream = 3 possibilities
        self.assertEqual(len(assemblies), 1 * 3 * 1)

        for assembly_idx, part in enumerate(promoter_collection):
            expected_sequence = ''.join([
                str(five_prime_flank.sequence),
                str(part.sequence),
                str(three_prime_flank.sequence)
            ])
            self.assertEqual(str(assemblies[assembly_idx].sequence), expected_sequence)
