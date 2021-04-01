"""Define the mechanism of SuperGSLFunction and an AST pass to invoke those functions."""
from typing import Optional, List
from inspect import getdoc

from supergsl.core.types import SuperGSLType
from supergsl.core.provider import SuperGSLProvider
from supergsl.core.symbol_table import SymbolTable
from supergsl.core.backend import DepthFirstNodeFilteredPass
from supergsl.core.exception import FunctionInvokeError, FunctionNotFoundError

#pylint: disable=E1136

class SuperGSLFunction(SuperGSLProvider, SuperGSLType):
    """Add a callable function to SuperGSL."""

    name: Optional[str] = None
    arguments : List[SuperGSLType] = []
    return_type : Optional[SuperGSLType] = None

    def resolve_import(self,
        symbol_table : SymbolTable,
        identifier : str,
        alias : str
    ) -> None:
        """Resolve the import of a function from this provider.

        """
        if identifier != self.name:
            raise FunctionNotFoundError('Function {} not provided by {}'.format(
                identifier, self))

        symbol_table.insert(alias or identifier, self)

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

    def execute(self, sgsl_args, child_nodes=None):
        """Called when the function is invoke in SuperGSL."""
        pass


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
            raise FunctionInvokeError('"%s" Return type does not match expectation. Expected: "%s", Actual: "%s"' % (
                node.function,
                expected_return_type,
                type(result_node[0])
            ))

        return result_node
