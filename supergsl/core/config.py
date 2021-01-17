import os
import logging
import json
from json import JSONDecodeError
from supergsl.core.exception import ConfigurationException

settings = {}
logger = logging.getLogger(__name__)

def load_config_file(config_file_path, optional=False):
    logger.debug('Attempting to load config file: "%s"..' % config_file_path)

    try:
        new_settings = json.load(open(config_file_path, 'r'))
        settings.update(new_settings)
    except FileNotFoundError as error:
        logger.debug('Failed to load config file: "%s".' % config_file_path)
        if not optional:
            raise ConfigurationException('Error loading config file: "%s"' % config_file_path, error)
    except JSONDecodeError as error:
        raise ConfigurationException('Error parsing config file: "%s"' % config_file_path, error)


base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/../'

load_config_file('%s/supergsl-config-default.json' % base_dir, optional=True)
load_config_file('%s/supergsl-config.json' % base_dir, optional=False)
