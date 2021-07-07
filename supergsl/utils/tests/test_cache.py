import pickle
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import Mock, patch
from pathlib import Path

from supergsl.utils.cache import FileCache

class FileCacheTestCase(TestCase):
    """Test the `FileCache` functionality."""

    def setUp(self):
        self.cache = FileCache('BOOM')

    @patch('supergsl.utils.cache.get_local_cache_path')
    def test_get_cached_path(self, get_local_cache_path_mock):
        """Local files are stored in the right place."""
        get_local_cache_path_mock.return_value = Path('/tmp','hi')
        result = self.cache.get_cached_path('HELLO')

        self.assertEqual(result, Path('/tmp', 'hi', 'HELLO.p'))

    def test_cache_miss_raises_key_error(self):
        """A item missing from the cache raises a `KeyError`"""
        temp_dir = TemporaryDirectory()
        temp_file_path = Path(temp_dir.name + '/test.p')

        self.cache.get_cached_path = Mock(return_value=temp_file_path)

        self.assertRaises(KeyError, self.cache.get, 'BOOOOMMM')

    def test_store_and_get_cached_part(self):
        temp_dir = TemporaryDirectory()
        temp_file_path = Path(temp_dir.name + '/test.p')

        self.cache.get_cached_path = Mock(return_value=temp_file_path)

        self.assertRaises(KeyError, self.cache.get, 'WHOOP')

        self.cache.store('WHOOP', {'HELLO': 'WORLD!'})

        result = self.cache.get('WHOOP')
        self.assertEqual(result, {
            'HELLO': 'WORLD!'
        })
        self.assertEqual(
            pickle.load(open(temp_file_path, 'rb')), {
                'HELLO': 'WORLD!'
            })
