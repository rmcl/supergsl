"""Implementation of part provider base class and plugin."""
from typing import Optional, List, Dict, Mapping, NamedTuple
from collections import namedtuple

from Bio.Seq import Seq

from supergsl.utils.resolve import import_class

from supergsl.core.plugin import SuperGSLPlugin
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.types import SuperGSLType
from supergsl.core.types.part import Part
from supergsl.core.types.position import Slice

from supergsl.core.exception import ConfigurationError, PartNotFoundError
from supergsl.core.constants import THREE_PRIME
from supergsl.core.provider import ProviderConfig


class PartProvider(SuperGSLProvider):

    def __init__(self, name : str, config : ProviderConfig):
        self._provider_name = name
        self.config = config

    def list_parts(self):
        """Return all parts available through this provider."""
        raise NotImplementedError(
            f'List parts is not supported by "{self.provider_name}" part provider.')

    def resolve_import(
        self,
        identifier : str,
        alias : Optional[str]
    ) -> Mapping[str, SuperGSLType]:
        """Resolve a part from the provider and register it in the symbol table."""
        part_identifier = alias or identifier

        # TODO: Maybe this is not solution, but temporary fix to get access to provider.
        if identifier == self._provider_name:
            return {
                part_identifier: self
            }

        return {
            part_identifier: self.get(identifier)
        }


    def get_child_part_by_slice(
        self,
        parent_part : Part,
        identifier : str,
        part_slice : Slice,
    ) -> Part:
        """Return a new part which is the child of the supplied parent."""
        raise NotImplementedError('Subclass to implement')


class ConstantPartDetail(NamedTuple):
    description : str
    sequence : str
    roles: List[str]


class ConstantPartProvider(PartProvider):
    """Load Parts from a python dictionary.

    This class can either be subclassed to provide constant sequences or
    instantiated with a `sequences` argument with desired parts.

    Provider Arguments:
    `name`: The name of the provider. This is used in your superGSL import statement.
    `provider_class`: For this provider, should be: "supergsl.core.parts.ConstantPartProvider".
    `sequences` (Mapping[str, ConstantPartDetail]): A dictionary containing the desired
        parts

    To utilize this class, pass the `sequence` argument or subclass and define
    a DEFAULT_PART_DETAILS dictionary.
    ```
    {
        <part-identifier>: ConstantPartDetail(
            <description>
            <sequence>
            [<list of roles>]
        )
    }
    ```
    """

    PART_DETAILS : Mapping[str, ConstantPartDetail] = {
        'part-name': ConstantPartDetail('part description', 'ATGC', ['LIST OF ROLES']),
    }

    def __init__(self, name : str, config : ProviderConfig):
        self._provider_name = name
        self._sequence_store = config.sequence_store
        self._config = config

        self.setup_part_details()
        self._cached_parts: Dict[str, Part] = {}

    def setup_part_details(self):
        if 'sequences' not in self._config.settings:
            self._part_details = self.PART_DETAILS
        else:
            self._part_details = self._config.settings['sequences']


    def list_parts(self):
        """Return all parts available through this provider."""
        return [
            self.get_part(identifier)
            for identifier in self._part_details.keys()
        ]

    def get_part_details(self, part_identifier):
        """Return constant details about a part."""
        return self._part_details[part_identifier]

    def get_default_part(self) -> Part:
        """Return the default part returned for `import <provider>` syntax.."""
        return self.get_part('default')

    def get(self, identifier : str) -> Part:
        """Return a part by identifier."""
        return self.get_part(identifier)

    def get_part(self, identifier : str) -> Part:
        """Retrieve a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        Return: `Part`
        """

        try:
            return self._cached_parts[identifier]
        except KeyError:
            pass

        try:
            description, reference_sequence, roles = \
                self.get_part_details(identifier)

        except KeyError as error:
            raise PartNotFoundError(
                'Part "%s" not found.' % identifier) from error

        if not isinstance(reference_sequence, Seq):
            reference_sequence = Seq(reference_sequence)

        sequence_entry = self._sequence_store.add_from_reference(reference_sequence)
        # TODO: Add annotations to this sequence., [
        #    Annotation(SEQ_ANNOTATION_ROLE, role)
        #    for role in roles
        #])

        part = Part(
            identifier,
            sequence_entry,
            provider=self,
            description=description,
            roles=roles)

        self._cached_parts[identifier] = part
        return part


class PartProviderPlugin(SuperGSLPlugin):
    """Plugin stub to register part providers."""
    name = 'part_provider'

    def register(self, compiler_settings):
        """Instantiate and register each part_providers defined in settings."""

        self.register_available_provider('constant_parts', ConstantPartProvider)

        if 'part_providers' not in compiler_settings:
            part_provider_configs = []
            print(
                'WARNING: '
                'No part providers have been defined. Check your supergGSL settings.')
        else:
            part_provider_configs = compiler_settings['part_providers']


        sequence_store = self.symbol_table.lookup('sequences')

        for provider_config in part_provider_configs:
            print('Initializing "%s"' % provider_config['name'])
            provider_class = import_class(provider_config['provider_class'])

            config = ProviderConfig(sequence_store, provider_config)
            provider_inst = provider_class(provider_config['name'], config)

            self.register_provider(provider_config['name'], provider_inst)
