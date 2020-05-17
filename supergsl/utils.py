import importlib


def import_class(class_path_str):
    """Import a class via str."""
    module_path, class_name = class_path_str.rsplit('.', 1)
    module = importlib.import_module(module_path)

    return getattr(module, class_name)
