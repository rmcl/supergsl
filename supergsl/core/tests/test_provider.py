"""Unit tests for the SuperGSL Provider base classes."""
from unittest import TestCase
from unittest.mock import Mock
from supergsl.core.provider import ProviderGroup


class ProviderGroupTestCase(TestCase):
    """Test that the parser correctly tokens into a valid AST."""
    maxDiff = None

    def test_provider_help(self):
        """Test that provider group help synthesizes help from all providers in group."""
        provider_group = ProviderGroup()

        for index in range(3):
            provider = Mock()
            provider.help = 'LET ME HELP YOU! %d' % index

            provider_group.add_provider(provider)

        self.assertEqual(
            provider_group.help, (
            'A collection of providers\n\n'
            'Mock\nLET ME HELP YOU! 0\n\n'
            'Mock\nLET ME HELP YOU! 1\n\n'
            'Mock\nLET ME HELP YOU! 2\n\n'))
