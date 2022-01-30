"""Test the Entrez part provider."""
import os
from unittest import TestCase
from unittest.mock import Mock, patch
from io import StringIO
from tempfile import TemporaryDirectory
from Bio import SeqIO
from supergsl.core.exception import ConfigurationError
from supergsl.plugins.ncbi.entrez import EntrezPartProvider


class EntrezPartProviderTestCase(TestCase):
    """Test case for the entrez part provider."""

    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        self.config = Mock()
        self.config.settings = {
            'local_genome_path': self.temp_dir.name + 'sgsl-lib/',
            'entrez_email': 'test@testexample.notreal',
            'efetch_args': {
                'id': 'boom123',
                'db': 'thedb'
            }
        }
        self.entrez_provider = EntrezPartProvider('test', self.config)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_check_for_required_settings(self):
        """ConfigurationError is raised if efetch_args and email not provided."""
        bad_config = Mock()
        bad_config.settings = {
            'efetch_args': {
                'id': 'boom123',
                'db': 'thedb'
            }
        }

        self.assertRaises(
            ConfigurationError,
            EntrezPartProvider,
            'test',
            bad_config)

    def test_get_local_file_path(self):
        """The local path should be generated and containing directory created."""
        full_path = self.entrez_provider.get_local_file_path()

        containing_dir = os.path.abspath(os.path.join(full_path, os.pardir))
        self.assertTrue(os.path.exists(containing_dir))

        self.assertEqual(
            full_path.replace(self.temp_dir.name, ''),
            'sgsl-lib/thedb-boom123.gb')

    @patch('supergsl.plugins.ncbi.entrez.Entrez')
    def test_retrieve_sequence_file(self, entrez_mock):
        test_gb_contents = 'TESTGBCONTENTS!'
        entrez_mock.efetch.return_value = StringIO(test_gb_contents)

        self.entrez_provider.retrieve_sequence_file()

        self.assertEqual(entrez_mock.email, self.config.settings['entrez_email'])
        entrez_mock.efetch.assert_called_once_with(
            id='boom123', db='thedb', rettype='gb', retmode='binary')
        self.assertEqual(
            open(self.entrez_provider.get_local_file_path()).read(),
            test_gb_contents)
