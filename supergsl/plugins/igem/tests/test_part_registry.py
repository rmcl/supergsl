"""Tests for the output core module."""
import requests
from unittest import TestCase
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.constants import SO_HOMOLOGOUS_REGION
from supergsl.plugins.igem.part_registry import (
    PartRegistry,
    BioBrickPartProvider
)
from supergsl.plugins.igem.types import BioBrickPart
from supergsl.plugins.igem.constants import (
    BIOBRICK_PREFIX_SEQUENCE,
    BIOBRICK_SUFFIX_SEQUENCE
)

class PartRegistryTestCase(TestCase):
    """Test the behavior of iGEM Part Registry provider"""
    maxDiff = None

    def setUp(self):
        self.mock_settings = {
            'repository_url': 'http://example.testbla/repo',
        }

    @patch('requests.get')
    def test_get_biobrick_part(self, request_mock_get):
        """Patch the http call to test the entire provider and part retrieval."""
        provider = PartRegistry('igem', self.mock_settings)

        request_mock_get.return_value.status_code = 200
        request_mock_get.return_value.content = open(
            'supergsl/plugins/builtin/providers/tests/BBa_E0040.xml', 'rb').read()

        part = provider.get_part('BBa_E0040')

        request_mock_get.assert_called_once_with(
            'http://example.testbla/repo/BBa_E0040/sbol',
            headers={
                'X-authorization': '',
                'Accept': 'text/plain'
            })

        self.assertEqual(type(part), BioBrickPart)
        self.assertEqual(part.identifier, 'BBa_E0040')
        self.assertEqual(part.description,
            'green fluorescent protein derived from jellyfish Aequeora victoria '
            'wild-type GFP (SwissProt: P42212')
        self.assertEqual(set(part.roles), set([
            'http://wiki.synbiohub.org/wiki/Terms/igem#partType/Coding',
            'http://identifiers.org/so/SO:0000316'
        ]))

        expected_sequence = ''.join([
            str(BIOBRICK_PREFIX_SEQUENCE).lower(),
            'atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgtta'
            'atgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaa'
            'atttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatggt'
            'gttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccg'
            'aaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagt'
            'caagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatgga'
            'aacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaac'
            'aaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagc'
            'agaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctg'
            'tccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttg'
            'taacagctgctgggattacacatggcatggatgaactatacaaataataa',
            str(BIOBRICK_SUFFIX_SEQUENCE).lower()
        ])

        self.assertEqual(
            str(part.sequence),
            expected_sequence)
        self.assertEqual(part.provider, provider)


class BioBrickPartProviderTestCase(TestCase):
    """Test that the BioBrickProvider can return prefix and suffix sequences."""

    def test_get_prefix_part(self):
        """Return a part with the biobrick prefix."""
        provider = BioBrickPartProvider('biobrick', {})

        prefix_part = provider.get_part('prefix')
        self.assertEqual(prefix_part.identifier, 'prefix')
        self.assertEqual(prefix_part.sequence, BIOBRICK_PREFIX_SEQUENCE)
        self.assertEqual(prefix_part.roles, [SO_HOMOLOGOUS_REGION])

    def test_get_suffix_part(self):
        """Return a part with the biobrick suffix."""
        provider = BioBrickPartProvider('biobrick', {})

        prefix_part = provider.get_part('suffix')
        self.assertEqual(prefix_part.identifier, 'suffix')
        self.assertEqual(prefix_part.sequence, BIOBRICK_SUFFIX_SEQUENCE)
        self.assertEqual(prefix_part.roles, [SO_HOMOLOGOUS_REGION])
