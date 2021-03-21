import importlib
import logging
import textwrap


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
