"""Support SuperGSL's symbol provider mechanism."""
from typing import List, Optional, Mapping
from supergsl.core.exception import NotFoundError
from supergsl.core.types import SuperGSLType
from supergsl.core.sequence import SequenceStore


class ProviderConfig:
    """Store config parameters for Part Providers"""

    def __init__(self, sequence_store : SequenceStore, settings : dict):
        self._sequence_store = sequence_store
        self._settings = settings

    @property
    def settings(self) -> dict:
        return self._settings

    @property
    def sequence_store(self) -> SequenceStore:
        return self._sequence_store


class SuperGSLProvider(SuperGSLType):
    """Base class to define objects that can provide symbols."""

    @property
    def provider_name(self):
        """Return the name of this part provider."""
        return self._provider_name

    @property
    def help(self):
        return ''

    def get(self, identifier : str) -> SuperGSLType:
        """Return a part by identifier.

        Arguments:
            identifier  A identifier to select a part from this provider
        """
        raise NotImplementedError('Subclass to implement')

    def save(self, part : SuperGSLType):
        """Save an object to the provider"""
        raise NotImplementedError('This provider does not support writing new parts.')

    def resolve_import(
        self,
        identifier : str,
        alias : Optional[str]
    ) -> Mapping[str, SuperGSLType]:
        """Resolve an identifier and return it"""
        raise NotImplementedError('Subclass to implement')


class ProviderGroup(SuperGSLProvider):
    """A group of providers allowing for sharing of the same module path."""

    _providers : List[SuperGSLProvider]

    def __init__(self):
        self._providers = []

    def add_provider(self, provider: SuperGSLProvider):
        """Add a provider to the group. Providers are searched in the order which they are added."""
        self._providers.append(provider)

    def get(self, identifier : str) -> SuperGSLType:
        """Return a part by identifier."""
        for provider in self._providers:
            try:
                return provider.get(identifier)
            except NotFoundError:
                continue
            except NotImplementedError:
                continue

        raise NotFoundError(f'{identifier} not found in module.')

    def __len__(self):
        """Return the number of providers in the group."""
        return len(self._providers)

    def __iter__(self):
        """Iterate through the providers in this group."""
        return iter(self._providers)

    def __contains__(self, provider : SuperGSLProvider):
        """Test if a provider is present in the group."""
        return provider in self._providers

    def __getitem__(self, index : int):
        """Return a provider by index."""
        return self._providers[index]

    @property
    def help(self) -> str:
        result = 'A collection of providers\n\n'
        for provider in self._providers:
            result += type(provider).__name__ + '\n'
            result += provider.help + '\n\n'
        return result

    def resolve_import(
        self,
        identifier : str,
        alias : Optional[str]
    ) -> Mapping[str, SuperGSLType]:
        """Import a identifier and register it in the symbol table."""

        for provider in self._providers:
            try:
                new_symbols = provider.resolve_import(identifier, alias)
            except NotFoundError:
                continue

            # Exit once we find a provider that doesn't raise a NotFoundError
            return new_symbols

        # If we try all the providers and don't get a result raise NotFoundError.
        raise NotFoundError('%s not found in module.' % (identifier))
