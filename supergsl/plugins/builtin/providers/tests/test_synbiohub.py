"""Tests for the output core module."""
import unittest
import requests
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.builtin.providers.synbiohub import SynBioHubPartProvider


class SynBioHubProviderTestCase(unittest.TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.mock_settings = {
            'repository_url': 'http://example.testbla/repo',
            'enable_part_cache' : False
        }

    def test_get_part_from_mocked_detail(self):
        """Confirm that the provider correctly initializes and returns a SuperGSL Part"""
        provider = SynBioHubPartProvider('igem', self.mock_settings)
        provider.get_part_details = Mock(
            return_value={
                'roles': [
                    'http://identifiers.org/so/SO:0000167',
                    'http://wiki.synbiohub.org/wiki/Terms/igem#partType/Regulatory',
                ],
                'description': 'constitutive promoter family member',
                'sequence': Seq('tttacggctagctcagtcctaggtatagtgctagc')
            })

        part = provider.get_part('BBa_J23106')

        provider.get_part_details.assert_called_once_with('BBa_J23106')
        self.assertEqual(part.identifier, 'BBa_J23106')
        self.assertEqual(part.description, 'constitutive promoter family member')
        self.assertEqual(part.roles, [
            'http://identifiers.org/so/SO:0000167',
            'http://wiki.synbiohub.org/wiki/Terms/igem#partType/Regulatory',
        ])
        self.assertEqual(str(part.sequence), 'tttacggctagctcagtcctaggtatagtgctagc')
        self.assertEqual(part.provider, provider)

    @patch('requests.get')
    def test_get_part(self, request_mock_get):
        """Patch the http call to test the entire provider and part retrieval."""
        provider = SynBioHubPartProvider('igem', self.mock_settings)

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

        self.assertEqual(part.identifier, 'BBa_E0040')
        self.assertEqual(part.description,
            'green fluorescent protein derived from jellyfish Aequeora victoria '
            'wild-type GFP (SwissProt: P42212')
        self.assertEqual(set(part.roles), set([
            'http://wiki.synbiohub.org/wiki/Terms/igem#partType/Coding',
            'http://identifiers.org/so/SO:0000316'
        ]))
        self.assertEqual(
            str(part.sequence),
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
            'taacagctgctgggattacacatggcatggatgaactatacaaataataa')
        self.assertEqual(part.provider, provider)
