"""Support SuperGSL's symbol provider mechanism."""
from typing import List
from supergsl.core.exception import NotFoundError
from supergsl.core.symbol_table import SymbolTable


class SuperGSLProvider:
    """Base class to define objects that can provide symbols."""

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Import a identifier and register it in the symbol table."""
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

    def resolve_import(
        self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Import a identifier and register it in the symbol table."""

        for provider in self._providers:
            try:
                provider.resolve_import(symbol_table, identifier, alias)
            except NotFoundError:
                continue

            # Exit once we find a provider that doesn't raise a NotFoundError
            return

        # If we try all the providers and don't get a result raise NotFoundError.
        raise NotFoundError('%s not found in module.' % (identifier))
