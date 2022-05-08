"""Define a minimal entrypoint to initialize the supergsl interpreter."""
from typing import Optional, List
from supergsl.core.config import load_settings
from .pipeline import CompilerPipeline


def get_interpreter(config_files : Optional[List[str]] = None):
    """Retrieve an instance of the SuperGSL compiler."""
    compiler_settings = load_settings(config_files)
    return CompilerPipeline(compiler_settings)
