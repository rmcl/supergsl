from unittest import TestCase
from unittest.mock import Mock
from Bio.Seq import Seq

from supergsl.core.exception import PartNotFoundError
from supergsl.core.parts.provider import PartProvider, ConstantPartProvider

class PartProviderTestCase(TestCase):
    """Test the `PartProvider` functionality."""

    def test_provider_name(self):
        """Provider name returns the expected name."""
        provider = PartProvider('BOOMPROVIDER', {})
        self.assertEqual(provider.provider_name, 'BOOMPROVIDER')

    def test_resolve_import_adds_parts_to_symbol_table(self):
        """Resolve import should add a part to the symbol table."""
        symbol_table = Mock()

        provider = PartProvider('BOOMPROVIDER', {})
        provider.get_part = Mock(return_value='PART!')

        new_symbols = provider.resolve_import('PARTIDENT', 'ALIASME')

        provider.get_part.assert_called_once_with('PARTIDENT')
        self.assertEqual(new_symbols['ALIASME'], 'PART!')

class ConstantPartProviderTestCase(TestCase):
    """Test case for `ConstantPartProvider`."""

    def test_constant_get_part(self):
        """Test that we can retrieve constant parts as `Part`."""

        provider = ConstantPartProvider('constant-part', {})

        provider.PART_DETAILS = {
            'test-part': (
                'THE GREATEST PART EVER',
                'ATGCAAATAGACAA',
                ['ROLE1', 'ROLE2']
            )
        }

        part = provider.get_part('test-part')

        self.assertEqual(part.identifier, 'test-part')
        self.assertEqual(part.sequence, 'ATGCAAATAGACAA')
        self.assertEqual(part.roles, ['ROLE1', 'ROLE2'])

    def test_constant_get_part_does_not_exist(self):
        """PartNotFoundError should be raised if part provider doesn't have the part."""
        provider = ConstantPartProvider('constant-part', {})

        with self.assertRaises(PartNotFoundError):
            provider.get_part('nonexistantpart')
