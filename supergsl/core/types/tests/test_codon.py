from unittest import TestCase
from supergsl.core.types.codon import CodonFrequencyTable


class CodonFrequencyTestCase(TestCase):
    """Test case for CodonFrequencyTable."""

    def test_print_codon_frequency_table(self):
        """Test that the repr output of a codon frequency table is as expected."""
        ex1 = {
            '*': {'TAA': 0.43, 'TAG': 0.18, 'TGA': 0.39},
            'A': {'GCA': 0.31, 'GCC': 0.2, 'GCG': 0.13, 'GCT': 0.36},
            'C': {'TGC': 0.45, 'TGT': 0.55},
        }

        table = CodonFrequencyTable('awesome', ex1).print()
        expected = '\n'.join([
            'Codon Frequencies for: awesome',
            'AmAcid\tCodon\tFrequency',
            '*\tTAA\t0.430000',
            '*\tTAG\t0.180000',
            '*\tTGA\t0.390000',
            '',
            'A\tGCA\t0.310000',
            'A\tGCC\t0.200000',
            'A\tGCG\t0.130000',
            'A\tGCT\t0.360000',
            '',
            'C\tTGC\t0.450000',
            'C\tTGT\t0.550000\n'
        ])
        self.assertEqual(str(table), expected)
