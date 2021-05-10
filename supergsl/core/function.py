"""Define the base class for SuperGSL Functions."""
from typing import Optional, List, Type, Dict, Tuple, Union
from inspect import getdoc

from supergsl.core.types import SuperGSLType

#pylint: disable=E1136


class SuperGSLFunction(SuperGSLType):
    """Add a callable function to SuperGSL."""

    name: Optional[str] = None
    compiler_settings : Optional[Dict] = None

    return_type : Optional[Type[SuperGSLType]] = None

    def __init__(self, compiler_settings : dict):
        self.settings = compiler_settings

    @classmethod
    def get_name(cls):
        if cls.name:
            return cls.name

        raise NotImplementedError('sGSL function definitions must specify a "name" in "%s".' % cls)

    @classmethod
    def get_help(cls) -> Optional[str]:
        return getdoc(cls)

    @classmethod
    def get_arguments(cls) -> List[Tuple[str, Type]]:
        """Return a list of expected arguments."""
        return []

    @classmethod
    def get_return_type(cls) -> Union[Type[SuperGSLType], Type[None]]:
        """Return the expected return value of the function."""
        return cls.return_type or type(None)

    def execute(self, params : dict):
        """Called when the function is invoke in SuperGSL."""
        raise NotImplementedError('Subclass to implement.')


class SuperGSLFunctionDeclaration(SuperGSLType):
    def __init__(self, function_class : Type[SuperGSLFunction], compiler_settings : dict):
        self.function_class = function_class
        self.compiler_settings = compiler_settings

    def eval(self, ast_node) -> SuperGSLFunction:
        return self.function_class(self.compiler_settings)
