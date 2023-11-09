"""Tests for the output core module."""
import unittest
import requests
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import SequenceStore
from supergsl.core.provider import ProviderConfig
from supergsl.plugins.builtin.providers.fasta import FastaPartProvider


class FASTAProviderTestCase(unittest.TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.temp_dir = TemporaryDirectory()
        temp_file_path = Path(self.temp_dir.name + '/test.fa')

        with open(temp_file_path, 'w+') as file_pointer:
            file_pointer.write('''>HELLO\natgc\n>YOYO\nttttaaaatgaacaaaa\n''')

        config = ProviderConfig(SequenceStore(), {
            'fasta_file_path': temp_file_path
        })
        self.provider = FastaPartProvider('bloop', config)

    def test_load_fasta_file_and_list_parts(self):
        """Load the temporary FASTA file and list all available parts."""
        self.assertEqual(
            self.provider.list_parts(),
            ['HELLO', 'YOYO'])

    def test_load_fasta_get_part(self):
        """Test retrieve a part from the fasta file and check its name and sequence."""
        part = self.provider.get_part('YOYO')
        self.assertEqual(part.identifier, 'YOYO')
        self.assertEqual(part.sequence, 'ttttaaaatgaacaaaa')
