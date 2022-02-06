"""Test the genbank file format part provider."""
from unittest import TestCase
from Bio import SeqIO
from supergsl.core.sequence import SequenceStore
from supergsl.core.provider import ProviderConfig
from supergsl.plugins.builtin.providers.genbank import GenBankFilePartProvider


class GenBankFilePartProviderTestCase(TestCase):
    """Test case for the file part provider."""

    def setUp(self):
        config = ProviderConfig(SequenceStore(), {
            'sequence_file_path': self.EXAMPLE_GB_FILE_PATH,
            'feature_types': ['CDS']
        })
        self.provider = GenBankFilePartProvider('test', config)

        self.expected_sequences = SeqIO.index(
            self.EXAMPLE_FASTA_FILE_PATH, 'fasta')


    def test_load_file_and_list_parts(self):
        """Test that we can load a genbank file and inspect the parts loaded."""
        parts = self.provider.list_parts()

        self.assertEqual(len(parts), 3)
        parts_by_identifier = {
            part.identifier: part
            for part in parts
        }

        capsid = parts_by_identifier['capsid']
        self.assertEqual(capsid.sequence, self.expected_sequences.get('capsid').seq)

        rep = parts_by_identifier['rep']
        self.assertEqual(rep.sequence, self.expected_sequences.get('rep').seq)

    EXAMPLE_GB_FILE_PATH = 'supergsl/plugins/ncbi/tests/fixtures/nucleotide-U89790.gb.gz'
    EXAMPLE_FASTA_FILE_PATH = 'supergsl/plugins/ncbi/tests/fixtures/nucleotide-U89790.fa'
