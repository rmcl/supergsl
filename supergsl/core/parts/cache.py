from typing import Any, Optional
import pickle
from datetime import datetime, timedelta
from pathlib import Path

from supergsl.utils import get_local_cache_path


class LocalFileCachePartProviderMixin(object):
    CACHED_FILE_EXPIRATION = timedelta(days=7)
    CACHED_FILE_EXTENSION = 'p'
    DELETE_EXPIRED_FILES = False

    _provider_name : str

    def get_cached_path(self, identifier : str) -> Path:
        """Return the path to the cache file."""
        filename = '%s.%s' % (identifier, self.CACHED_FILE_EXTENSION)
        return Path(
            get_local_cache_path(self._provider_name),
            filename)

    def cached_file_exists(self, cached_file_path : Path) -> bool:
        """Return true if a cached file exists and has not expired."""
        if not cached_file_path.exists():
            return False

        last_modified_time = datetime.fromtimestamp(cached_file_path.stat().st_mtime)
        if datetime.now() - last_modified_time > self.CACHED_FILE_EXPIRATION:
            if self.DELETE_EXPIRED_FILES:
                # if the file has expired remove it to avoid confusion next time.
                cached_file_path.unlink()
            return False

        return True

    def store_part_details(self, identifier : str, data : Any):
        """Store part details in the cache."""
        cached_file_path = self.get_cached_path(identifier)
        with open(cached_file_path, 'wb+') as file_handle:
            pickle.dump(data, file_handle)

    def get_cached_part_details(self, identifier) -> Any:
        """Retrieve part details from the cache."""
        cached_file_path = self.get_cached_path(identifier)
        if self.cached_file_exists(cached_file_path):
            with open(cached_file_path, 'rb') as file_handle:
                return pickle.load(file_handle)
        else:
            result = self.get_part_details(identifier)
            self.store_part_details(identifier, result)
            return result

    def get_part_details(self, identifier : str) -> Any:
        """Retrieve part details from the source."""
        raise NotImplementedError('Subclass to implement')
