import importlib
import logging
import textwrap
from pathlib import Path
from supergsl.core.config import load_settings

def import_class(class_path_str):
    """Import a class via str."""
    module_path, class_name = class_path_str.rsplit('.', 1)
    module = importlib.import_module(module_path)

    return getattr(module, class_name)


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
