import os
import json
import importlib
import logging
import textwrap
from pathlib import Path


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

def setup_supergsl_config_files():
    conf_dir = Path('~/.supergsl')
    conf_dir = conf_dir.expanduser()
    lib_dir = conf_dir.joinpath('sgsl-lib/')
    conf_file_path = conf_dir.joinpath('config.json')
    example_config_file = Path(__file__).parent.joinpath('../supergsl-config.json.example')

    try:
        Path.mkdir(conf_dir)
    except FileExistsError:
        print(
            'Cannot setup config files at "%s" because that directory '
            'already exists.' % conf_dir)
        return

    Path.mkdir(lib_dir, exist_ok=True)

    json_config = json.load(open(example_config_file))
    json_config['lib_dir'] = lib_dir.as_posix()

    with open(conf_file_path, 'w+') as conf_file_fp:
        json.dump(json_config, conf_file_fp, indent=4)
        conf_file_fp.write('\n')
