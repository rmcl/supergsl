"""Test the genbank file format part provider."""
from unittest import TestCase
from Bio import SeqIO
from supergsl.plugins.ncbi.genbank_provider import GenBankFilePartProvider


class GenBankFilePartProviderTestCase(TestCase):
    """Test case for the file part provider."""

    def setUp(self):
        self.provider = GenBankFilePartProvider('test', {
            'sequence_file_path': self.EXAMPLE_GB_FILE_PATH,
            'feature_types': ['CDS']
        })

        self.expected_sequences = SeqIO.index(
            self.EXAMPLE_FASTA_FILE_PATH, 'fasta')


    def test_load_file_and_list_parts(self):
        """Test that we can load a genbank file and inspect the parts loaded."""
        parts = self.provider.list_parts()

        self.assertEqual(len(parts), 2)
        parts_by_identifier = {
            part.identifier: part
            for part in parts
        }

        capsid = parts_by_identifier['capsid']
        self.assertEqual(capsid.start.get_absolute_position_in_reference()[1], 2259)
        self.assertEqual(capsid.end.get_absolute_position_in_reference()[1], 4464)
        self.assertEqual(capsid.get_sequence().seq, self.expected_sequences.get('capsid').seq)

        rep = parts_by_identifier['rep']
        self.assertEqual(rep.start.get_absolute_position_in_reference()[1], 371)
        self.assertEqual(rep.end.get_absolute_position_in_reference()[1], 2243)
        self.assertEqual(rep.get_sequence().seq, self.expected_sequences.get('rep').seq)

    EXAMPLE_GB_FILE_PATH = 'supergsl/plugins/ncbi/test/fixtures/nucleotide-U89790.gb.gz'
    EXAMPLE_FASTA_FILE_PATH = 'supergsl/plugins/ncbi/test/fixtures/nucleotide-U89790.fa'
