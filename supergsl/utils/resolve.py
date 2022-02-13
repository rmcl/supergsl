import importlib
from typing import Optional, List, Union, Callable
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.exception import ConfigurationError, ProviderNotFoundError


def import_class(class_path_str : str):
    """Import a class via str."""
    module_path, class_name = class_path_str.rsplit('.', 1)
    module = importlib.import_module(module_path)

    return getattr(module, class_name)

def resolve_provider_import(
    provider : Union[str, Callable],
    symbol_table : SymbolTable
):
    """Determine if the provider is a class or a string that needs the be imported."""
    # If provider is callable then we are done.
    if callable(provider):
        return provider

    # If provider is not callable then attempt to import it.
    try:
        return import_class(provider)
    except ValueError:
        pass

    # If we can't determine the provider class then maybe it is already been
    # registered. Attempt to find it in the symbol table.
    available_providers = symbol_table.enter_nested_scope('available_imports')
    try:
        return available_providers[provider]
    except KeyError as error:
        raise ProviderNotFoundError(f'Unknown Provider "{provider}"') from error


def resolve_import(
    symbol_table : SymbolTable,
    module_path : List[str],
    identifier : str,
    alias : Optional[str]
):
    """Resolve an import at a particular module path in the given symbol table."""
    import_table = symbol_table.enter_nested_scope('imports')

    module_path_str = '.'.join(module_path)
    provider = import_table.lookup(module_path_str)
    if not isinstance(provider, SuperGSLProvider):
        raise ConfigurationError('"%s" is not a provider. It is a %s' % (
            module_path_str,
            type(provider)
        ))

    new_symbols = provider.resolve_import(
        identifier,
        alias)

    for symbol_identifier, symbol_value in new_symbols.items():
        symbol_table.insert(symbol_identifier, symbol_value)

    return new_symbols
