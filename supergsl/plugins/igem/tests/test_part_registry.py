"""Tests for the output core module."""
import requests
from unittest import TestCase
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.core.sequence import SequenceStore
from supergsl.core.provider import ProviderConfig
from supergsl.core.types.role import SO_HOMOLOGOUS_REGION
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
        sequence_store = SequenceStore()
        self.mock_config = ProviderConfig(sequence_store, {
            'repository_url': 'http://example.testbla/repo',
        })

class BioBrickPartProviderTestCase(TestCase):
    """Test that the BioBrickProvider can return prefix and suffix sequences."""

    def setUp(self):
        sequence_store = SequenceStore()
        self.mock_config = ProviderConfig(sequence_store, {
            'repository_url': 'http://example.testbla/repo',
        })

    def test_get_prefix_part(self):
        """Return a part with the biobrick prefix."""
        provider = BioBrickPartProvider('biobrick', self.mock_config)

        prefix_part = provider.get_part('prefix')
        self.assertEqual(prefix_part.identifier, 'prefix')
        self.assertEqual(prefix_part.sequence, BIOBRICK_PREFIX_SEQUENCE)
        self.assertEqual(prefix_part.roles, [SO_HOMOLOGOUS_REGION])

    def test_get_suffix_part(self):
        """Return a part with the biobrick suffix."""
        provider = BioBrickPartProvider('biobrick', self.mock_config)

        prefix_part = provider.get_part('suffix')
        self.assertEqual(prefix_part.identifier, 'suffix')
        self.assertEqual(prefix_part.sequence, BIOBRICK_SUFFIX_SEQUENCE)
        self.assertEqual(prefix_part.roles, [SO_HOMOLOGOUS_REGION])
