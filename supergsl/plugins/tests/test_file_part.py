from unittest import TestCase
from unittest.mock import Mock
from Bio import SeqIO
from supergsl.core.pipeline import CompilerPipeline
from supergsl.tests.fixtures import SuperGSLIntegrationFixtures
from supergsl.plugins.file_part import FeatureTableWithFastaPartProvider
from supergsl.core.parts.provider import PartProviderConfig
from supergsl.core.sequence import SequenceStore


class FilePartTestCases(TestCase):

    def setUp(self):
        self.maxDiff = None
        self.expected_sequences = SeqIO.index(
            'supergsl/tests/expected_sequences.fasta', 'fasta')

        self.store = SequenceStore()

        self.fixtures = SuperGSLIntegrationFixtures()
        settings = self.fixtures.get_supergsl_settings()
        config = PartProviderConfig(self.store, settings['part_providers'][0])

        self.provider = FeatureTableWithFastaPartProvider('TESTER', config)

    def test_get_gene_forward_strand(self):
        """Test that the provider properly returns a gene on the forward strand of a genome."""
        self.provider.load()

        entry, gene_feature = self.provider.get_gene('GAL3')
        del gene_feature['Notes']
        self.assertEqual(gene_feature, {
            '': '',
            '?': '',
            'aliases': 'GAL3,transcriptional regulator GAL3',
            'chrom#': '4',
            'chromname': 'chromosome 4',
            'from': 463433,
            'gene': 'GAL3',
            'id': 'YDR009W',
            'qualifier': '?',
            'strand': 'W',
            'systematic': 'YDR009W',
            'to': 464996,
            'type': 'gene'
        })

        self.assertEqual(entry.sequence, self.expected_sequences['gGAL3'].seq)


    def test_get_gene_reverse_strand(self):
        """Test that the provider properly returns a gene on the reverse strand of a genome."""
        self.provider.load()

        entry, gene_feature = self.provider.get_gene('HO')
        del gene_feature['Notes']
        self.assertEqual(gene_feature, {
            '': '',
            '?': '',
            'aliases': 'HO',
            'chrom#': '4',
            'chromname': 'chromosome 4',
            'from': 48031,
            'gene': 'HO',
            'id': 'YDL227C',
            'qualifier': '?',
            'strand': 'C',
            'systematic': 'YDL227C',
            'to': 46270,
            'type': 'gene'
        })

        self.assertEqual(entry.sequence, self.expected_sequences['gHO'].seq)
