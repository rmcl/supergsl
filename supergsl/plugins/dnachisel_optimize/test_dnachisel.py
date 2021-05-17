from unittest import TestCase
from mock import Mock
from Bio import SeqIO
from Bio.Seq import Seq
from supergsl.core.types.builtin import AminoAcidSequence
from supergsl.plugins.dnachisel_optimize import DNAChiselOptimizeFunction

from dnachisel import reverse_translate

class DnaChiselTestCases(TestCase):

    def setUp(self):
        self.maxDiff = None

        settings = {}
        self.optimize = DNAChiselOptimizeFunction(settings)

        """
        self.expected_sequences = SeqIO.index(
            'supergsl/test/expected_sequences.fasta', 'fasta')

        self.fixtures = SuperGSLIntegrationFixtures()
        self.compiler_settings = self.fixtures.get_supergsl_settings(
            extra_plugins=[
                'supergsl.plugins.dnachisel_optimize'
            ]
        )
        """

    def test_execute(self):
        """Call dnachisel codon_optimize execute function."""
        params = {
            'aa_sequence': AminoAcidSequence(Seq('MAAATCAGAGAAAAC')),
            'num_results': 2
        }
        self.optimize.create_new_sequence = Mock(
            return_value=Seq('ATGAAAC'))
        results = self.optimize.execute(params)

        result_sequences = [
            item.get_sequence()
            for item in results
        ]

        self.assertEqual(result_sequences, [
            Seq('ATGAAAC'),
            Seq('ATGAAAC')
        ])

    def test_create_new_sequence(self):
        """Run the dnachisel optimizer and get a new DNA sequence."""
        target_protein = 'MAAATCAGAGAAAAC'
        naive_target_sequence = reverse_translate(target_protein)
        result = self.optimize.create_new_sequence(
            naive_target_sequence,
            None,
            []
        )

        self.assertEqual(
            Seq(result).translate(),
            target_protein)
