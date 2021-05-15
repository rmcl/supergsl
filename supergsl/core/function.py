"""Define the base class for SuperGSL Functions."""
from typing import (
    Optional,
    List,
    Type,
    Dict,
    Tuple,
    Union,
    Any
)
from inspect import getdoc

from supergsl.core.types import SuperGSLType
from supergsl.core.exception import FunctionInvokeError
#pylint: disable=E1136


class SuperGSLFunction(SuperGSLType):
    """Add a callable function to SuperGSL."""

    name: Optional[str] = None
    compiler_settings : Optional[Dict] = None

    arguments : List[Tuple[str, Type]] = []
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
        return cls.arguments

    @classmethod
    def get_return_type(cls) -> Union[Type[SuperGSLType], Type[None]]:
        """Return the expected return value of the function."""
        return cls.return_type or type(None)

    def execute(self, params : dict):
        """Called when the function is invoke in SuperGSL."""
        raise NotImplementedError('Subclass to implement.')

    ### Helper Methods for Function Invocation

    def check_function_result(self, result) -> None:
        """Check the result of a function against the declared return type."""
        expected_return_type = self.get_return_type()
        if not expected_return_type:
            expected_return_type = type(None)

        if not isinstance(result, expected_return_type):
            raise FunctionInvokeError(
                '"%s" Return type does not match expectation. Expected: "%s", Actual: "%s"' % (
                    self,
                    expected_return_type,
                    type(result)
                ))

    def build_argument_map(
        self,
        positional_arguments,
        child_arguments
    ) -> Dict[str, Any]:
        """Build an argument dict from positional arugments."""

        expected_arguments = self.get_arguments()
        if len(positional_arguments) != len(expected_arguments):
            raise FunctionInvokeError(
                'Number of positional arguments does not match function definition. '
                'Expected %d, but received %d' % (
                    len(expected_arguments),
                    len(positional_arguments))
            )

        function_parameters = {}
        for argument_idx, argument_details in enumerate(expected_arguments):
            argument_key, expected_argument_type = argument_details
            argument_value = positional_arguments[argument_idx]

            if argument_key == 'children':
                raise Exception('Cannot define an argument named "children". It is reserved.')

            if not isinstance(argument_value, expected_argument_type):
                raise FunctionInvokeError(
                    'Provided type does not match expectation. '
                    'Expected %s, but received %s' % (
                        expected_argument_type,
                        type(argument_value))
                )

            function_parameters[argument_key] = argument_value

        return function_parameters

    def evaluate_arguments_and_execute(
        self,
        positional_arguments : List[Any],
        child_arguments : List[Any]
    ):
        """Perform type checking and execute the desired SuperGSL Function."""

        arguments = self.build_argument_map(
            positional_arguments,
            child_arguments)

        result = self.execute(arguments)
        self.check_function_result(result)

        return result




class SuperGSLFunctionDeclaration(SuperGSLType):
    def __init__(self, function_class : Type[SuperGSLFunction], compiler_settings : dict):
        self.function_class = function_class
        self.compiler_settings = compiler_settings

    def eval(self) -> SuperGSLFunction:
        return self.function_class(self.compiler_settings)
