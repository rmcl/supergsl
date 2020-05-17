import json
from supergsl.core.exception import ConfigurationException

settings = {}
try:
    settings = json.load(open('supergsl-config.json', 'r'))
except Exception as error:
    raise ConfigurationException('Error loading config', error)
