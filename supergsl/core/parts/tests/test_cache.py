import pickle
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import Mock, patch
from pathlib import Path
from Bio.Seq import Seq

from supergsl.core.exception import PartNotFoundError
from supergsl.core.parts.cache import LocalFileCachePartProviderMixin

class LocalFileCachePartProviderMixinTestCase(TestCase):
    """Test the `LocalFileCachePartProviderMixin` functionality."""

    def setUp(self):
        self.provider = LocalFileCachePartProviderMixin()
        self.provider._provider_name = 'BOOM'

    @patch('supergsl.core.parts.cache.get_local_cache_path')
    def test_get_cached_path(self, get_local_cache_path_mock):
        get_local_cache_path_mock.return_value = Path('/tmp','hi')
        result = self.provider.get_cached_path('HELLO')

        self.assertEqual(result, Path('/tmp', 'hi', 'HELLO.p'))

    @patch('supergsl.core.parts.cache.get_local_cache_path')
    def test_get_cached_part_details(self, get_local_cache_path_mock):
        temp_dir = TemporaryDirectory()
        temp_file_path = Path(temp_dir.name + '/test.p')

        self.provider.get_cached_path = Mock(return_value=temp_file_path)
        self.provider.get_part_details = Mock(return_value={
            'HELLO': 'WORLD!'
        })

        result = self.provider.get_cached_part_details('WHOOP')
        self.assertEquals(result, {
            'HELLO': 'WORLD!'
        })
        self.assertEqual(
            pickle.load(open(temp_file_path, 'rb')), {
                'HELLO': 'WORLD!'
            })

        self.provider.get_cached_part_details('WHOOP')
        self.provider.get_cached_part_details('WHOOP')

        self.provider.get_part_details.assert_called_once_with('WHOOP')
