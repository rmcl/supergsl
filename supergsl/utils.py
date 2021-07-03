import importlib
import logging
import textwrap
from pathlib import Path
from typing import cast, Optional, List
from supergsl.core.config import load_settings
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.provider import SuperGSLProvider

def import_class(class_path_str : str):
    """Import a class via str."""
    module_path, class_name = class_path_str.rsplit('.', 1)
    module = importlib.import_module(module_path)

    return getattr(module, class_name)


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
        raise Exception('"%s" is not a provider. It is a %s' % (
            module_path_str,
            type(provider)
        ))

    provider.resolve_import(
        symbol_table,
        identifier,
        alias)


def get_logger(class_inst):
    full_path = ".".join([
        class_inst.__class__.__module__,
        class_inst.__class__.__name__
    ])

    return logging.getLogger(full_path)

def get_local_cache_path(provider_name) -> Path:
    settings = load_settings()
    cache_path = Path(
        settings['local_cache_path'],
        provider_name)

    cache_path.mkdir(parents=True, exist_ok=True)

    return cache_path

def get_logo():
    """Return the SuperGSL ASCII logo.

    >>> from supergsl.utils import get_logo
    >>> print(get_logo())

    """
    return textwrap.dedent(
        """
          ____                         ____ ____  _
         / ___| _   _ _ __   ___ _ __ / ___/ ___|| |
         \___ \| | | | '_ \ / _ \ '__| |  _\___ \| |
          ___) | |_| | |_) |  __/ |  | |_| |___) | |___
         |____/ \__,_| .__/ \___|_|   \____|____/|_____|
                     |_|

        """)


def display_symbol_table(symbol_table, depth=0):
    """Recursively display the table and all of its nested scopes."""
    indent = ' ' * (depth*4)
    print()
    print(indent + 'Table: {} ----'.format(symbol_table.name))
    for key, val in symbol_table._symbols.items():
        print(indent + '    Key: {}, Type: {}'.format(
            key,
            type(val)
        ))

    for nested_scope in symbol_table._nested_scopes.values():
        display_symbol_table(nested_scope, depth+1)

    print(indent + '--- End Table: %s ----' % symbol_table.name)
