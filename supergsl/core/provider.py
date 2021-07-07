"""Support SuperGSL's symbol provider mechanism."""
from typing import List, Optional, Mapping
from supergsl.core.exception import NotFoundError
from supergsl.core.types import SuperGSLType


class SuperGSLProvider:
    """Base class to define objects that can provide symbols."""

    @property
    def help(self):
        return ''

    def resolve_import(
        self,
        identifier : str,
        alias : Optional[str]
    ) -> Mapping[str, SuperGSLType]:
        """Resolve an identifier and return it"""
        raise NotImplementedError('Subclass to implement')


class ProviderGroup(SuperGSLProvider):
    """A group of providers allowing for sharing of the same module path."""

    def __init__(self):
        self._providers : List[SuperGSLProvider] = []

    def add_provider(self, provider: SuperGSLProvider):
        """Add a provider to the group. Providers are searched in the order which they are added."""
        self._providers.append(provider)

    def __contains__(self, provider : SuperGSLProvider):
        """Test if a provider is present in the group."""
        return provider in self._providers

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
