import importlib
import logging


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
