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
from typeguard import check_type
from inspect import getdoc

from supergsl.core.types import SuperGSLType
from supergsl.core.exception import FunctionInvokeError
from supergsl.core.sequence import SequenceStore
from supergsl.core.provider import ProviderConfig

#pylint: disable=E1136


class SuperGSLFunction(SuperGSLType):
    """Add a callable function to SuperGSL."""

    name: Optional[str] = None
    compiler_settings : Optional[Dict] = None

    arguments : List[Tuple[str, Type]] = []
    return_type : Optional[Type[SuperGSLType]] = None

    def __init__(self, config : ProviderConfig):
        self.settings = config.settings
        self.sequence_store = config.sequence_store

    @classmethod
    def get_name(cls):
        if cls.name:
            return cls.name

        raise NotImplementedError('sGSL function definitions must specify a "name" in "%s".' % cls)

    @property
    def help(cls) -> Optional[str]:
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

        try:
            check_type(self.name, result, expected_return_type)
        except TypeError as error:
            raise FunctionInvokeError(
                '"%s" Return type does not match expectation. Expected: "%s", Actual: "%s"' % (
                    self,
                    expected_return_type,
                    type(result)
                )) from error

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
                raise FunctionInvokeError('Cannot define an argument named "children". It is reserved.')

            try:
                check_type(self.name, argument_value, expected_argument_type)
            except TypeError as error:
                raise FunctionInvokeError(
                    'Provided type does not match expectation. '
                    'Expected %s, but received %s' % (
                        expected_argument_type,
                        type(argument_value))
                ) from error

            function_parameters[argument_key] = argument_value

        if child_arguments:
            function_parameters['children'] = child_arguments

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
    """Store details related to the declaration of a SuperGSL Function."""

    def __init__(self, function_class : Type[SuperGSLFunction], compiler_settings : dict):
        self.function_class = function_class
        self.compiler_settings = compiler_settings
        self.sequence_store : Optional[SequenceStore] = None

    def set_sequence_store(self, sequence_store : SequenceStore):
        self.sequence_store = sequence_store

    def eval(self) -> SuperGSLFunction:
        if not self.sequence_store:
            raise Exception('SequenceStore not set.')

        config = ProviderConfig(self.sequence_store, self.compiler_settings)
        return self.function_class(config)
