from re import Pattern
from typing import Optional, Tuple, List, Dict, Type
from collections import OrderedDict
from supergsl.core.exception import SymbolNotFoundError, ProviderNotFoundError, NotFoundError
from supergsl.core.provider import SuperGSLProvider


class SymbolTable(object):

    def __init__(self):
        self._import_path_to_providers : Dict[str, List[Type[SuperGSLProvider]]] = {}
        self._symbols_providers: List[Tuple[Pattern, SuperGSLProvider]] = []
        self._symbols = {}

    def register(self, provider_import_path: str, provider_class: Type[SuperGSLProvider]) -> None:
        """Register a provider to be available for import a given import path."""
        try:
            self._import_path_to_providers[provider_import_path].append(provider_class)
        except KeyError:
            self._import_path_to_providers[provider_import_path] = [provider_class]

    def get_providers_for_path(self, provider_import_path) -> List[Type[SuperGSLProvider]]:
        """Return all providers associated with an import path."""
        providers = self._import_path_to_providers.get(provider_import_path, None)
        if not providers:
            raise ProviderNotFoundError('Provider for %s not found.' % provider_import_path)

        return providers

    def resolve_symbol(
        self,
        import_path: str,
        import_name: str,
        alias: Optional[str] = None
    ):
        """Resolve a symbol from an import statement of a superGSL program."""

        print('Attempting to resolve symbol: %s, %s, %s' % (import_path, import_name, alias))

        active_alias = alias or import_name
        try:
            self.get_symbol(active_alias)
        except SymbolNotFoundError:
            pass
        else:
            raise Exception('Alias "%s" imported twice.' % active_alias)

        provider_classes = self.get_providers_for_path(import_path)
        for provider_class in provider_classes:
            try:
                alias_pattern = provider_class.resolve_import(import_name, active_alias)
            except NotFoundError:
                continue

            self._symbols_providers.append((alias_pattern, provider_class))
            return self.get_symbol(active_alias)

        raise NotFoundError('Symbol "{}" could not be resolved in path "{}"'.format(
            import_name,
            import_path
        ))

    def get_symbol(self, identifier):
        """Retrieve an imported symbol."""

        try:
            return self._symbols[identifier]
        except KeyError:
            pass

        for symbol_pattern, symbol_handler in self._symbols_providers:
            match = symbol_pattern.fullmatch(identifier)
            if match:
                self._symbols[identifier] = symbol_handler.get_symbol(match)
                return self._symbols[identifier]

        raise SymbolNotFoundError('"%s" is not defined.' % identifier)
