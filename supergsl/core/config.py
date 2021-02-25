import logging
import json
from json import JSONDecodeError
from supergsl.core.exception import ConfigurationException

logger = logging.getLogger(__name__)


def load_config_file(config_file_path, optional=False):
    logger.debug('Attempting to load config file: "%s"..' % config_file_path)

    try:
        return json.load(open(config_file_path, 'r'))
    except FileNotFoundError as error:
        logger.debug('Failed to load config file: "%s".' % config_file_path)
        if not optional:
            raise ConfigurationException('Error loading config file: "%s"' % config_file_path, error)
    except JSONDecodeError as error:
        raise ConfigurationException('Error parsing config file: "%s"' % config_file_path, error)

def load_settings():
    """Load typical SuperGSL Config files."""
    settings = load_config_file('supergsl-config-default.json', optional=True) or {}
    settings.update(load_config_file('supergsl-config.json', optional=False))

    return settings
