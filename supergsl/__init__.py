from typing import Optional, List
from .core.config import load_settings
from .core.pipeline import CompilerPipeline


def get_compiler(config_files : Optional[List[str]] = None):
    """Retrieve an instance of the SuperGSL compiler."""
    compiler_settings = load_settings(config_files)
    return CompilerPipeline(compiler_settings)
