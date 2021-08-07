from unittest import TestCase
from Bio.Seq import Seq

from supergsl.utils.primers import build_two_part_assembly_primers
from supergsl.core.types.primer import Primer


class PrimerUtilTestCases(TestCase):
    """Test the primer utility methods."""

    def test_build_two_part_assembly_primers(self):
        """Assembly primers should create expected primers and sequence."""

        five_prime_part_seq = Seq('atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg')
        three_primer_part_seq = Seq('ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc')
        payload_seq = Seq('ATGAACAT')

        result = build_two_part_assembly_primers(
            five_prime_part_seq, three_primer_part_seq, payload_seq)

        self.assertEqual(result['assembly_sequence'], Seq(''.join([
            str(five_prime_part_seq),
            str(payload_seq),
            str(three_primer_part_seq)
        ])))

        expected_primers = {
            'terminal_forward_primer': 'gaagggttagcagtcat',
            'terminal_reverse_primer': 'ggcaacttgaccaag',
            'internal_forward_primer': 'ctggtgggtttggATGTTCATcatcgtaagtttcg',
            'internal_reverse_primer': 'cctggtgggtttggATGTTCATcatcgtaagtttcgaacga',
        }
        for primer_name, expected_primer_seq in expected_primers.items():
            self.assertEqual(result[primer_name].sequence, Seq(expected_primer_seq),
                '%s does not match expected sequence' % primer_name)
