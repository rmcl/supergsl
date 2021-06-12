"""Tests for the SBOL output module."""
import unittest
import requests
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from supergsl.core.tests.fixtures import SuperGSLCoreFixtures
from supergsl.plugins.synbiohub.provider import SynBioHubPartProvider


class SBOLOutputTestCase(unittest.TestCase):
    """Test the behavior of SynBioHub provider"""
    maxDiff = None

    def setUp(self):
        self.mock_settings = {
            'repository_url': 'http://example.testbla/repo',
        }

    def test_get_part_from_mocked_detail(self):
        pass
