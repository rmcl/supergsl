"""Define the mechanism of SuperGSLFunction and an AST pass to invoke those functions."""
from typing import Optional, List, Type, Dict
from inspect import getdoc

from supergsl.core.types import SuperGSLType
from supergsl.core.provider import SuperGSLProvider
#from supergsl.core.exception import FunctionInvokeError, FunctionNotFoundError

#pylint: disable=E1136


class SuperGSLFunction(SuperGSLType):
    """Add a callable function to SuperGSL."""

    name: Optional[str] = None
    compiler_settings : Optional[Dict] = None

    return_type : Optional[SuperGSLType] = None

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
    def get_arguments(cls):
        """Return a list of expected arguments."""
        return cls.arguments

    @classmethod
    def get_return_type(cls):
        """Return the expected return value of the function."""
        return cls.return_type

    def execute(self, params : dict):
        """Called when the function is invoke in SuperGSL."""
        raise NotImplementedError('Subclass to implement.')


class SuperGSLFunctionDeclaration(SuperGSLProvider):
    def __init__(self, function_class : Type[SuperGSLFunction], compiler_settings : dict):
        self.function_class = function_class
        self.compiler_settings = compiler_settings

    def eval(self, ast_node) -> SuperGSLFunction:
        return self.function_class(self.compiler_settings)

'''
class InvokeFunctionPass(DepthFirstNodeFilteredPass):
    """Traverse the AST and execute encountered SuperGSLFunctions."""

    def get_node_handlers(self):
        return {
            'FunctionInvocation': self.visit_function_invoke_node,
        }

    def visit_function_invoke_node(self, node):
        print('INVOKE', node.params, node.identifier)

        if node.params is not None:
            print('WARNING!!! PASSING PARAMS NOT IMPLEMENTED YET!!!!!!!')
            #raise NotImplementedError('Passing parameters to functions is not yet implemented.')

        args = {}
        result_node = node.function.execute(args, node.get_definition_list())
        expected_return_type = node.function.get_return_type()
        if not isinstance(result_node, expected_return_type):
            raise FunctionInvokeError(
                '"%s" Return type does not match expectation. Expected: "%s", Actual: "%s"' % (
                    node.function,
                    expected_return_type,
                    type(result_node[0])
                ))

        return result_node
'''
