import os
import logging
import json
from typing import Optional, List
from json import JSONDecodeError
from supergsl.core.exception import ConfigurationError, ConfigFileNotFound

logger = logging.getLogger(__name__)

# Global module variable to cache loaded settings.
__cached_settings = None

def load_config_file(config_file_path, optional=False):
    logger.debug('Attempting to load config file: "%s"..' % config_file_path)

    try:
        return json.load(open(config_file_path, 'r'))
    except FileNotFoundError as error:
        logger.debug('Failed to load config file: "%s".' % config_file_path)
        if not optional:
            raise ConfigFileNotFound(f'Config file not found: "{config_file_path}"')
    except JSONDecodeError as error:
        raise ConfigurationError('Error parsing config file: "%s"' % config_file_path, error)

def load_settings(config_files : Optional[List[str]] = None):
    """Load typical SuperGSL Config files.

    Supply an optional list of paths to config files to override SuperGSL default.
    """
    global __cached_settings
    if __cached_settings is not None:
        return __cached_settings

    if config_files:
        print('Loading config files: %s' % config_files)
        settings = {}
        for config_file in config_files:
            settings.update(load_config_file(config_file, optional=False))

    else:
        settings = load_config_file('supergsl-config-default.json', optional=True) or {}
        try:
            settings.update(load_config_file('supergsl-config.json', optional=False) or {})
        except ConfigFileNotFound:
            print(
                'WARNING: supergsl-config.json not found. Using default config options. '
                'Very limited providers and plugin functionality will be available.')

    if 'local_cache_path' not in settings:
        settings['local_cache_path'] = './sgsl-lib/'

    __cached_settings = settings
    return settings
